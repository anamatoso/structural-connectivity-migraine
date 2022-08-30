#!/bin/bash
set -e
# Have to be in data directory
MAINDIR="/strombolihome/amatoso/tese/data"
display_usage() {
	echo "$(basename $0) [Subject and type of session]"
	# echo "This script uses FSL to analyze diffusion data. It requires X arguments: 
	# 	1) Subject DWI Directory (string)- has to have data.nii.gz (...b0b1b2000.nii.gz), nodif_brain_mask.nii.gz, bvecs (...b0b1000b2000), bvals
	#	example: "sub-control019"
	#	2) Type of session. Example: ses-midcycle
	}

	if [ $# -le 0 ] # if there are 0 argument or less
	then
		display_usage
		exit 1
	fi

SUBDIR="/strombolihome/amatoso/tese/data/${1}" #example name: sub-control019_ses-midcycle
cd $SUBDIR
if [ ! -d "bedpostdata.bedpostX" ] 
then 
	########################### STEP 1 ###################################
	#             		  Prepare for bedpostx                           #
	######################################################################

	# create bedpost directory inside subject and copy files to there

	mkdir -p bedpostdata
	cp *_${SUBDIR}_clean.nii.gz ./bedpostdata/data.nii.gz
	cp *_${SUBDIR}_clean_mask.nii.gz ./bedpostdata/nodif_brain_mask.nii.gz
	cp "/home/amatoso/tese/data/bvals_132dir.bval" ./bedpostdata/bvals
	cp *_${SUBDIR}_eddy.eddy_rotated_bvecs.bvec ./bedpostdata/bvecs


	# check if its all good. It should be.
	# bedpostx_datacheck bedpostdata

	########################### STEP 1 ###################################
	#             Fiber orientation estimation                           #
	######################################################################

	# Creates the directory needed for the probablistic tractography
	bedpostx ./bedpostdata -model 2 -n 4
	cd $MAINDIR
else 
	cd $MAINDIR
	exit 0
fi