#!/bin/bash
set -e
# Have to be in data directory
MAINDIR=$(pwd)
display_usage() {
	echo "$(basename $0) [Subject and type of session]"
	# echo "This script uses FSL to analyze diffusion data. It requires X arguments: 
	# 	1) Subject DWI Directory (string)- has to have data.nii.gz (...b0b1b2000.nii.gz), nodif_brain_mask.nii.gz, bvecs (...b0b1000b2000), bvals
	#	example: "sub-control019"
	#	2) Type of session. Example: ses-midcycle
	}

	if [ $# -le 0 ] # if there are 1 argument or less
	then
		display_usage
		exit 1
	fi
MAINDIR2="/strombolihome/amatoso/tese/data"
DIR=$1 #example name: sub-control019_ses-midcycle

########################### STEP 1 ###################################
#             Prepare for atlas for probtrackx                       #
######################################################################

SUBDIR="${MAINDIR2}/${DIR}"
cd $SUBDIR

# divide the atlas into several files each one with one ROI
/home/amatoso/tese/divide_atlas_seedgmwm.sh "${SUBDIR}/mrtrix_outputs_bvals2/gmwmseed_atlas.mif" $SUBDIR
ROIS="${MAINDIR2}/${DIR}/list_rois_gmwminterface.txt"


mrconvert "${MAINDIR}/${DIR}/mrtrix_outputs_bvals2/gmwmseed_atlas.mif" "${MAINDIR}/${DIR}/mrtrix_outputs_bvals2/gmwmseed_atlas.nii.gz" -force

########################### STEP 3 ###################################
#                        Tractography                                #
######################################################################

n_fibers=5000 		# number of fibers
fiber_len=250 		# length of fiber in mm
step=1 				# step of each guess
mask="${SUBDIR}/bedpostdata.bedpostX/nodif_brain_mask.nii.gz"
seed="${SUBDIR}/bedpostdata.bedpostX/merged"
output_dir="${SUBDIR}/probtrackX_outputs_omat3"
gmwm_mask="${MAINDIR2}/${DIR}/mrtrix_outputs_bvals2/gmwmseed_atlas.nii.gz"


probtrackx2 --network -x $ROIS --usef=0.05 -l --onewaycondition --omatrix3 --target3=$ROIS -c 0.2 -S $fiber_len --steplength=$step --stop=$gmwm_mask -P $n_fibers --fibthresh=0.01 --distthresh=4 --sampvox=0.0 --forcedir --opd -s $seed -m $mask --dir=$output_dir

# copy connectivity matrix to folder
cd $output_dir
cp fdt_network_matrix "${MAINDIR2}/matrix_data/${DIR}_fsl_matrix_bval2_omat3_seedgmwminterface"
cd ..
rm -rf $output_dir
# go back to data directory
cd $MAINDIR
