#!/bin/bash

set -e # Stop on errors

display_usage() {
	echo "This script uses FSL's probtrackx2 to create a connectivity matrix."
	echo "Usage: $(basename $0) [Subject and type of session]"
	echo "t requires 1 argument: the subject directory. Example: sub-control019_ses-midcycle"
	}

	if [ $# -le 0 ] # if there are 1 argument or less
	then
		display_usage
		exit 1
	fi

# Have to be in data directory
MAINDIR=$(pwd)

DIR=$1 #example name: sub-control019_ses-midcycle

########################### STEP 1 ###################################
#             Prepare for atlas for probtrackx                       #
######################################################################

SUBDIR="${MAINDIR}/${DIR}"
cd $SUBDIR

# divide the atlas into several files each one with one ROI
#/home/amatoso/tese/divide_atlas.sh "${SUBDIR}/mrtrix_outputs_bvals2/atlas.mif" $SUBDIR -> done in MRtrix pipeline

ROIS="${MAINDIR}/${DIR}/list_rois.txt"

# Convert GM/WM boundary to nii.gz if it isn't already
if [ ! -f "${MAINDIR}/${DIR}/mrtrix_outputs_bvals2/gmwmseed_atlas.nii.gz" ]; then
	mrconvert "${MAINDIR}/${DIR}/mrtrix_outputs_bvals2/gmwmseed_atlas.mif" "${MAINDIR}/${DIR}/mrtrix_outputs_bvals2/gmwmseed_atlas.nii.gz"
fi

########################### STEP 3 ###################################
#                        Tractography                                #
######################################################################
n_fibers=5000 		# number of fibers
fiber_len=250 		# length of fiber in mm
step=1 				# step of each guess
mask="${SUBDIR}/bedpostdata.bedpostX/nodif_brain_mask.nii.gz"
seed="${SUBDIR}/bedpostdata.bedpostX/merged"
output_dir="${SUBDIR}/probtrackX_outputs_omat3"
gmwm_mask="${MAINDIR}/${DIR}/mrtrix_outputs_bvals2/gmwmseed_atlas.nii.gz"

probtrackx2 --network -x $ROIS --usef=0.05 -l --onewaycondition --omatrix3 --target3=$ROIS -c 0.2 -S $fiber_len --steplength=$step --stop=$gmwm_mask -P $n_fibers --fibthresh=0.01 --distthresh=4 --sampvox=0.0 --forcedir --opd -s $seed -m $mask --dir=$output_dir

# copy connectivity matrix to matrices folder
cp "${output_dir}/fdt_network_matrix" "${MAINDIR}/matrix_data/${DIR}_fsl_matrix_bval2_omat3"

# delete output directory because it is memory heavy
rm -rf $output_dir

# go back to data directory
cd $MAINDIR
