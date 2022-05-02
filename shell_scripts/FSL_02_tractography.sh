#!/bin/bash

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

DIR=$1 #example name: sub-control019_ses-midcycle

########################### STEP 1 ###################################
#             Prepare for atlas for probtrackx                       #
######################################################################

SUBDIR="${MAINDIR}/${DIR}"
cd $SUBDIR

# divide the atlas into several files each one with one ROI
/home/amatoso/tese/divide_atlas.sh "${SUBDIR}/mrtrix_outputs/atlas.mif" $SUBDIR
ROIS="${SUBDIR}/list_rois.txt"

########################### STEP 3 ###################################
#                        Tractography                                #
######################################################################

n_fibers=5000 		# number of fibers
fiber_len=250 		# length of fiber in mm
step=1 				# step of each guess
mask="${MAINDIR}/${SUBDIR}/bedpostdata.bedpostX/nodif_brain_mask"
seed="${MAINDIR}/${SUBDIR}/bedpostdata.bedpostX/merged"
output_dir="${MAINDIR}/${SUBDIR}/probtrackX_outputs"

probtrackx2 --network -x $ROIS -l --onewaycondition --omatrix1 -c 0.2 -S $fiber_len --steplength=$step -P $n_fibers --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd -s $seed -m $mask --dir=$output_dir
 #probtrackx -s "basename" -m bet_mask -x gm_boundary --omatrix3 connmatrix --target3 mask_roi
 # probtrackx2 -x /home/amatoso/tese/data/sub-control019_ses-midcycle_altered/mrtrix_outputs_og/gmwmSeed_coreg_bin.nii.gz -l --onewaycondition --network --omatrix1 --omatrix3 --target3 /home/amatoso/tese/data/sub-control019_ses-midcycle/mrtrix_outputs/atlas.nii.gz -c 0.2 -S 250 --steplength=1 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd -s /home/amatoso/tese/data/sub-control019_ses-midcycle_altered/bedpostdata.bedpostX/merged -m /home/amatoso/tese/data/sub-control019_ses-midcycle_altered/bedpostdata.bedpostX/nodif_brain_mask --dir=/home/amatoso/tese/data/sub-control019_ses-midcycle_altered/probtrackX_outputs5000 

# change name of connectivity matrix
cd $output_dir
mv "${DIR}_fsl_matrix"

# go back to data directory
cd $MAINDIR
