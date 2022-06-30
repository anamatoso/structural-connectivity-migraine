#!/bin/bash
set -e
# Have to be in data directory
MAINDIR=$(pwd)
display_usage() {
        echo "$(basename $0) [Subject and type of session]"
        # echo "This script uses MRTrix to analyze diffusion data. It requires X arguments: 
        #       1) Subject DWI Directory. Example: sub-control019
        #       2) Type of session. Example: ses-midcycle
        }

        if [ $# -le 0 ] # if there are 1 argument or less
        then
                display_usage
                exit 1
        fi

DIR=$1 #example name: sub-control019_ses-midcycle
SUB=${DIR:0:14} # example name: sub-control019

########################### STEP 1 ###################################
#                         Prepare data and directories                                   #
######################################################################
ATLAS="/home/amatoso/tese/data/AAL116_intersect.nii.gz"
#SPINEMASK="/home/amatoso/tese/data/spinalcord_mask.nii.gz"
SUBDIR="${MAINDIR}/${DIR}" 
ANATDIR="${MAINDIR}/${SUB}" #example name: sub-control019
ANAT="${ANATDIR}/${SUB}_restored-MPRAGE_brain.nii.gz"
cd $SUBDIR

# WHEN DOING FOR REAL, CHECK IF THERE IS MORE DATA WITH MORE BVALUES
DWI="${SUBDIR}/*clean.nii.gz"
DWIMASK="${SUBDIR}/*clean_mask.nii.gz"
BVEC="${SUBDIR}/*rotated_bvecs.bvec"
BVAL="/home/amatoso/tese/data/bvals_132dir.bval"

mkdir -p mrtrix_outputs_bvals2
cd mrtrix_outputs_bvals2

########################### STEP 5 ###################################
#                         Coregister atlas to the data                   #
######################################################################

#Coregister atlas to struct space and convert to mrtrix format
applywarp -i $ATLAS -r $ANAT --out=atlasintersect_2struct --warp="${ANATDIR}/reg_nonlinear_invwarp_T1tostandard_2mm.nii.gz"
mrconvert atlasintersect_2struct.nii.gz atlasintersect_2struct.mif -force

#applywarp -i $SPINEMASK -r $ANAT --out=spinemask_2struct --warp="${ANATDIR}/reg_nonlinear_invwarp_T1tostandard_2mm.nii.gz"
#mrconvert spinemask_2struct.nii.gz spinemask_2struct.mif -force

#Coregister atlas to diff space
mrtransform atlasintersect_2struct.mif -linear diff2struct_mrtrix.txt -inverse atlasintersect_coreg.mif -force
#mrtransform spinemask_2struct.mif -linear diff2struct_mrtrix.txt -inverse spinemask_coreg.mif -force

# Make sure the values of the atlas are integer
mrcalc atlasintersect_coreg.mif -round -datatype uint32 atlasintersect.mif -force

# Restric the GM/WM boundary to the atlas
mrcalc atlasintersect.mif gmwmSeed_coreg.mif -mult gmwmseed_atlasintersect.mif -force

# Create file with sizes of ROIs (in voxels) - for use in matlab
mrconvert atlasintersect.mif atlasintersect.nii.gz -force
fslstats atlasintersect.nii.gz -h 117 >> roi_size_w0.txt
tail -n +2 roi_size_w0.txt >> roi_size_intersect.txt
rm roi_size_w0.txt


########################### STEP 6 ###################################
#                        Additional Steps                                     #
######################################################################

mrconvert $DWI dwi.mif -fslgrad $BVEC $BVAL -force
mrconvert $ANAT anat.mif -force


# Extract all five tissue catagories (1=GM; 2=Subcortical GM; 3=WM; 4=CSF; 5=Pathological tissue)
5ttgen fsl anat.mif 5tt_nocoreg.mif -premasked -nocrop -force

# The following series of commands will take the average of the b0 images (which have the best contrast), convert them and the 5tt image to NIFTI format, and use it for coregistration.
dwiextract dwi.mif - -bzero | mrmath - mean mean_b0_processed.mif -axis 3 -force
mrconvert mean_b0_processed.mif mean_b0_processed.nii.gz -force
mrconvert 5tt_nocoreg.mif 5tt_nocoreg.nii.gz -force

# Uses FSL commands fslroi and flirt to create a transformation matrix for regisitration between the tissue map and the b0 images
fslroi 5tt_nocoreg.nii.gz 5tt_vol0.nii.gz 0 1 #Extract the first volume of the 5tt dataset (since flirt can only use 3D images, not 4D images)

flirt -in mean_b0_processed.nii.gz -ref 5tt_vol0.nii.gz -interp nearestneighbour -dof 6 -omat diff2struct_fsl.mat
transformconvert diff2struct_fsl.mat mean_b0_processed.nii.gz 5tt_nocoreg.nii.gz flirt_import diff2struct_mrtrix.txt -force
mrtransform 5tt_nocoreg.mif -linear diff2struct_mrtrix.txt -inverse 5tt_coreg.mif -force

#######


########################### STEP 6 ###################################
#                 Run the streamline analysis                        #
######################################################################

# Create streamlines: maxlength=250mm, 10M seeds
tckgen -act 5tt_coreg.mif -algorithm SD_STREAM -seed_gmwmi gmwmseed_atlasintersect.mif -maxlength 250 -select 10000000 wmfod_norm.mif tracksintersect.tck -force

# Reduce the number of streamlines with tcksift
tcksift2 -act 5tt_coreg.mif tracksintersect.tck wmfod_norm.mif siftintersect.txt -force

########################### STEP 7 ###################################
#                         Creating the connectome                            #
######################################################################

#Creating the connectome 
tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in siftintersect.txt tracksintersect.tck atlasintersect.mif "${MAINDIR}/matrix_data/${DIR}_mrtrix_matrix_bval2_intersect.csv" -force

rm -f 5tt_nocoreg.mif 5tt_nocoreg.nii.gz 5tt_nocoreg.mif dwi.mif anat.mif mean_b0_processed.mif mean_b0_processed.nii.gz 5tt_vol0.nii.gz
cd $MAINDIR


