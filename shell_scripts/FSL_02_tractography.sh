#!/bin/bash

set -e #Stop on errors

display_usage() {
	tput bold 
	echo ""
	echo "Description"
	tput sgr0
	echo ""
	echo "This script uses FSL's probtrackx2 to create a connectivity matrix. Note that either the MRtrix pipeline needs to be done first or the command for the division of the atlas should be uncommented."
	echo "It uses the option omatrix1 in FSL"
	tput bold 
	echo ""
	echo "Usage"
	tput sgr0
	echo ""
	echo "./$(basename $0) [Subject and type of session]"
	echo "It requires 1 argument: the subject directory. Example: sub-control019_ses-midcycle"
	}

	if [ $# -le 0 ] # if there no argument
	then
		display_usage
		exit 1
	fi

# Have to be in data directory
MAINDIR=$(pwd)

DIR=$1

########################### STEP 1 ###################################
#             Prepare for atlas for probtrackx                       #
######################################################################

SUBDIR="${MAINDIR}/data/${DIR}"
cd $SUBDIR

# divide the atlas into several files each one with one ROI
# ${MAINDIR}/divide_atlas.sh "${SUBDIR}/mrtrix_outputs_bvals2/atlas.mif" $SUBDIR -> done in MRtrix pipeline

ROIS="${MAINDIR}/data/${DIR}/list_rois.txt"

# Convert GM/WM boundary to nii.gz if it isn't already
if [ ! -f "${MAINDIR}/data/${DIR}/mrtrix_outputs_bvals2/gmwmseed_atlas.nii.gz" ]; then
	mrconvert "${MAINDIR}/data/${DIR}/mrtrix_outputs_bvals2/gmwmseed_atlas.mif" "${MAINDIR}/data/${DIR}/mrtrix_outputs_bvals2/gmwmseed_atlas.nii.gz"
fi

########################### STEP 3 ###################################
#                        Tractography                                #
######################################################################
n_fibers=5000 		# number of fibers
fiber_len=250 		# length of fiber in mm
step=1 				# step of each guess in mm
mask="${SUBDIR}/bedpostdata.bedpostX/nodif_brain_mask.nii.gz"
seed="${SUBDIR}/bedpostdata.bedpostX/merged"
output_dir="${SUBDIR}/probtrackX_outputs"
gmwm_mask="${MAINDIR}/data/${DIR}/mrtrix_outputs_bvals2/gmwmseed_atlas.nii.gz"

probtrackx2 --network -x $ROIS -l --onewaycondition --omatrix1 -c 0.2 -S $fiber_len --steplength=$step --stop=$gmwm_mask -P $n_fibers --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd -s $seed -m $mask --dir=$output_dir

# Copy connectivity matrix to matrices folder
cp "${output_dir}/fdt_network_matrix" "${MAINDIR}/matrix_data/${DIR}_fsl_matrix_bval2"

# Delete output directory because it is memory heavy
rm -rf $output_dir

# Go back to data directory
cd $MAINDIR
