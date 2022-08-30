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

MAINDIR2="/strombolihome/amatoso/tese/data"
DIR=$1 #example name: sub-control019_ses-midcycle
SUB=${DIR:0:14} # example name: sub-control019

########################### STEP 1 ###################################
#                         Prepare data and directories                                   #
######################################################################
ATLAS="/strombolihome/amatoso/desikan.nii.gz"
SUBDIR="${MAINDIR2}/${DIR}"
ANATDIR="${MAINDIR2}/${SUB}" #example name: sub-control019
ANAT="${ANATDIR}/${SUB}_restored-MPRAGE_brain.nii.gz"
cd $SUBDIR

# WHEN DOING FOR REAL, CHECK IF THERE IS MORE DATA WITH MORE BVALUES
DWI="${SUBDIR}/*clean.nii.gz"
DWIMASK="${SUBDIR}/*clean_mask.nii.gz"
BVEC="${SUBDIR}/*rotated_bvecs.bvec"
BVAL="/strombolihome/amatoso/tese/data/bvals_132dir.bval"

mkdir -p mrtrix_outputs_bvals2
cd mrtrix_outputs_bvals2

########################### STEP 5 ###################################
#                         Coregister atlas to the data                   #
######################################################################

#Coregister atlas to struct space and convert to mrtrix format
applywarp -i $ATLAS -r $ANAT --out=atlas_2struct_desikan --warp="${ANATDIR}/reg_nonlinear_invwarp_T1tostandard_2mm.nii.gz"
mrconvert atlas_2struct_desikan.nii.gz atlas_2struct_desikan.mif -force

#Coregister atlas to diff space
mrtransform atlas_2struct_desikan.mif -linear diff2struct_mrtrix.txt -inverse atlas_coreg_desikan.mif -force
#mrtransform spinemask_2struct.mif -linear diff2struct_mrtrix.txt -inverse spinemask_coreg.mif -force

# Make sure the values of the atlas are integer
mrcalc atlas_coreg_desikan.mif -round -datatype uint32 atlas_desikan.mif -force

# Restric the GM/WM boundary to the atlas
mrcalc atlas_desikan.mif gmwmSeed_coreg.mif -mult gmwmseed_atlas_desikan.mif -force

# Create file with sizes of ROIs (in voxels) - for use in matlab
mrconvert atlas_desikan.mif atlas_desikan.nii.gz -force
fslstats atlasintersec_desikant.nii.gz -h 117 >> roi_size_w0.txt
tail -n +2 roi_size_w0.txt >> roi_size_desikan.txt
rm roi_size_w0.txt

########################### STEP 6 ###################################
#                 Run the streamline analysis                        #
######################################################################

# Create streamlines: maxlength=250mm, 10M seeds
tckgen -act 5tt_coreg.mif -algorithm SD_STREAM -seed_gmwmi gmwmseed_atlas_desikan.mif -maxlength 250 -select 10000000 wmfod_norm.mif tracks_desikan.tck -force

# Reduce the number of streamlines with tcksift
tcksift2 -act 5tt_coreg.mif tracks_desikan.tck wmfod_norm.mif sift_desikan.txt -force

########################### STEP 7 ###################################
#                         Creating the connectome                            #
######################################################################

#Creating the connectome 
tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in sift_desikan.txt tracks_desikan.tck atlas_desikan.mif "${MAINDIR2}/matrix_data/${DIR}_mrtrix_matrix_bval2_desikan.csv" -force

rm -f 5tt_nocoreg.nii.gz 5tt_nocoreg.mif mean_b0_processed.mif mean_b0_processed.nii.gz 5tt_vol0.nii.gz anat.mif atlas_2struct_desikan.nii.gz atlas_2struct_desikan.mif atlas_coreg_desikan.mif atlas_desikan.nii.gz atlas_coreg.mif gmfod.mif gmfod_norm.mif csffod.mif csffod_norm.mif mask.mif wmfod.mif dwi.mif wm.txt gm.txt csf.txt

cd $MAINDIR

