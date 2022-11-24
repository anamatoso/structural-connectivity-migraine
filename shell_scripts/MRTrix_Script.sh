#!/bin/bash

set -e # Stop on errors

display_usage() {
	echo "This script uses MRTrix to analyze diffusion data and create a connectivity matrix."
	echo "Usage: $(basename $0) [Directory]"
	echo "It requires 1 argument: Subject DWI Directory. Example: sub-control019_ses-midcycle"
	}

	if [ $# -le 0 ] # if there are 0
	then
		# can't run function
		display_usage
		exit 1
	fi

MAINDIR=$(pwd)

DIR=$1 #example name: sub-control019_ses-midcycle
SUB=${DIR:0:14} # crop name. result=sub-control019

########################### STEP 1 ###################################
#             		  Prepare data and directories					 #
######################################################################

# Get directories
SUBDIR="${MAINDIR}/${DIR}" 
ANATDIR="${MAINDIR}/${SUB}"
ANAT="${ANATDIR}/${SUB}_restored-MPRAGE_brain.nii.gz"
cd $SUBDIR

# Get data
ATLAS="AAL116_intersect.nii.gz"
DWI="${SUBDIR}/*clean.nii.gz"
DWIMASK="${SUBDIR}/*clean_mask.nii.gz"
BVEC="${SUBDIR}/*rotated_bvecs.bvec"
BVAL="bvals_132dir.bval"

# Make output folder and change to that
mkdir -p mrtrix_outputs_bvals2
cd mrtrix_outputs_bvals2

########################### STEP 2 ###################################
#	      		  Convert data to .mif format				   	     #
######################################################################

# Convert data do mif fomat
mrconvert $DWI dwi.mif -fslgrad $BVEC $BVAL -force
mrconvert $ANAT anat.mif -force

# Convert mask do mif fomat
mrconvert $DWIMASK mask.mif -force

########################### STEP 3 ###################################
#             Basis function for each tissue type                    #
######################################################################

# Create a basis function from the subject's DWI data. The "dhollander" function is best used for multi-shell acquisitions; it will estimate different basis functions for each tissue type. For single-shell acquisition, use the "tournier" function instead
dwi2response dhollander dwi.mif wm.txt gm.txt csf.txt -force

# Performs multishell-multitissue constrained spherical deconvolution, using the basis functions estimated above
dwi2fod msmt_csd dwi.mif -mask mask.mif wm.txt wmfod.mif gm.txt gmfod.mif csf.txt csffod.mif -force

# Creates an image of the fiber orientation densities overlaid onto the estimated tissues (Blue=WM; Green=GM; Red=CSF)
# You should see FOD's mostly within the white matter. These can be viewed later with the command "mrview vf.mif -odf.load_sh wmfod.mif"
#mrconvert -coord 3 0 wmfod.mif - | mrcat csffod.mif gmfod.mif - vf.mif

# Normalize the FODs to enable comparison between subjects
mtnormalise wmfod.mif wmfod_norm.mif gmfod.mif gmfod_norm.mif csffod.mif csffod_norm.mif -mask mask.mif -force

########################### STEP 4 ###################################
#            Create a GM/WM boundary for seed analysis               #
######################################################################

# Extract all five tissue catagories (1=GM; 2=Subcortical GM; 3=WM; 4=CSF; 5=Pathological tissue)
5ttgen fsl anat.mif 5tt_nocoreg.mif -premasked -nocrop -force #premasked because the data is already just the brain and no crop to maintain the size of the image

# The following series of commands will take the average of the b0 images (which have the best contrast), convert them and the 5tt image to NIFTI format, and use it for coregistration.
dwiextract dwi.mif - -bzero | mrmath - mean mean_b0_processed.mif -axis 3 -force
mrconvert mean_b0_processed.mif mean_b0_processed.nii.gz -force
mrconvert 5tt_nocoreg.mif 5tt_nocoreg.nii.gz -force

# Uses FSL commands fslroi and flirt to create a transformation matrix for regisitration between the tissue map and the b0 images
fslroi 5tt_nocoreg.nii.gz 5tt_vol0.nii.gz 0 1 #Extract the first volume of the 5tt dataset (since flirt can only use 3D images, not 4D images)

flirt -in mean_b0_processed.nii.gz -ref 5tt_vol0.nii.gz -interp nearestneighbour -dof 6 -omat diff2struct_fsl.mat
transformconvert diff2struct_fsl.mat mean_b0_processed.nii.gz 5tt_nocoreg.nii.gz flirt_import diff2struct_mrtrix.txt -force
mrtransform 5tt_nocoreg.mif -linear diff2struct_mrtrix.txt -inverse 5tt_coreg.mif -force

#Create a seed region along the GM/WM boundary
5tt2gmwmi 5tt_coreg.mif gmwmSeed_coreg.mif -force

########################### STEP 5 ###################################
#             		  Coregister atlas to the data                   #
######################################################################

#Coregister atlas to struct space and convert to mrtrix format
applywarp -i $ATLAS -r $ANAT --out=atlas_2struct --warp="${ANATDIR}/reg_nonlinear_invwarp_T1tostandard_2mm.nii.gz"
mrconvert atlas_2struct.nii.gz atlas_2struct.mif -force

#Coregister atlas to diffusion space
mrtransform atlas_2struct.mif -linear diff2struct_mrtrix.txt -inverse atlas_coreg.mif -force

# Make sure the values of the atlas are integer
mrcalc atlas_coreg.mif -round -datatype uint32 atlas.mif -force

# Restric the GM/WM boundary to the atlas
mrcalc atlas.mif gmwmSeed_coreg.mif -mult gmwmseed_atlas.mif -force

# for probtrackx (in FSL pipeline)
mrconvert gmwmseed_atlas.mif gmwmseed_atlas.nii.gz
./divide_atlas.sh "${SUBDIR}/mrtrix_outputs_bvals2/atlas.mif" $SUBDIR

# Create file with sizes of ROIs (in voxels) - to use in matlab to normalize matrices for ROI size
mrconvert atlas.mif atlas.nii.gz -force
fslstats atlas.nii.gz -h 117 >> roi_size_w0.txt
tail -n +2 roi_size_w0.txt >> roi_size.txt
rm roi_size_w0.txt

########################### STEP 6 ###################################
#                 Run the streamline generation                      #
######################################################################

# Create streamlines: maxlength=250mm, 10M seeds
tckgen -act 5tt_coreg.mif -seed_gmwmi gmwmseed_atlas.mif -maxlength 250 -select 10000000 wmfod_norm.mif tracks.tck -force

# Filter streamlines with tcksift
tcksift2 -act 5tt_coreg.mif -out_mu sift_mu.txt -out_coeffs sift_coeffs.txt tracks.tck wmfod_norm.mif sift.txt -force

########################### STEP 7 ###################################
#             		  Creating the connectome	                     #
######################################################################

#Creating the connectome 
tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in sift.txt tracks.tck atlas.mif "${MAINDIR}/matrix_data/${DIR}_mrtrix_matrix_bval2_intersect.csv" -force

# Remove big files and only keep the ones useful for the critical steps
rm -f tracks.tck 5tt_nocoreg.nii.gz 5tt_nocoreg.mif mean_b0_processed.mif mean_b0_processed.nii.gz 5tt_vol0.nii.gz anat.mif atlas_2struct_desikan.nii.gz atlas_2struct_desikan.mif atlas_coreg_desikan.mif atlas_desikan.nii.gz atlas_coreg.mif gmfod.mif gmfod_norm.mif csffod.mif csffod_norm.mif mask.mif wmfod.mif dwi.mif wm.txt gm.txt csf.txt

# Go back to main directory
cd $MAINDIR
