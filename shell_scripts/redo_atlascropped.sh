
set -e
# Have to be in data directory
MAINDIR=$(pwd)
display_usage() {
	echo "$(basename $0) [Subject and type of session]"
	# echo "This script uses MRTrix to analyze diffusion data. It requires X arguments: 
	# 	1) Subject DWI Directory. Example: sub-control019
	#	2) Type of session. Example: ses-midcycle
	}

	if [ $# -le 0 ] # if there are 1 argument or less
	then
		display_usage
		exit 1
	fi

DIR=$1 #example name: sub-control019_ses-midcycle

########################### STEP 1 ###################################
#             		  Prepare data and directories					 #
######################################################################

SUBDIR="${MAINDIR}/${DIR}" 
cd $SUBDIR
cd mrtrix_outputs_bvals2
#Coregister atlas to diff space
mrtransform atlas_2struct.mif -linear diff2struct_mrtrix.txt -inverse -template dwi.mif atlas_coreg.mif -force

# Make sure the values of the atlas are integer
mrcalc atlas_coreg.mif -round -datatype uint32 atlas.mif -force

# Restric the GM/WM boundary to the atlas
# mrcalc atlas.mif gmwmSeed_coreg.mif -mult gmwmseed_atlas.mif -force

# Create file with sizes of ROIs (in voxels) - for use in matlab
mrconvert atlas.mif atlas.nii.gz -force
fslstats atlas.nii.gz -h 117 >> roi_size_w0.txt
tail -n +2 roi_size_w0.txt >> roi_size.txt
rm roi_size_w0.txt

########################### STEP 6 ###################################
#                 Run the streamline analysis                        #
######################################################################

# Create streamlines: maxlength=250mm, 10M seeds
tckgen -act 5tt_coreg.mif -algorithm SD_STREAM -seed_gmwmi gmwmseed_atlas.mif -maxlength 250 -select 10000000 wmfod_norm.mif tracks.tck -force

# Reduce the number of streamlines with tcksift
tcksift2 -act 5tt_coreg.mif -out_mu sift_mu.txt -out_coeffs sift_coeffs.txt tracks.tck wmfod_norm.mif sift.txt -force

########################### STEP 7 ###################################
#             		  Creating the connectome	                     #
######################################################################

#Creating the connectome 
tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in sift.txt tracks.tck atlas.mif "${MAINDIR}/matrix_data/${DIR}_mrtrix_matrix_bval2.csv" -force

cd $MAINDIR
