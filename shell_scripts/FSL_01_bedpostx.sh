#!/bin/bash

# Have to be in data directory
MAINDIR=$(pwd) #data
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

SUBDIR=$1 #example name: sub-control019_ses-midcycle

########################### STEP 1 ###################################
#             		  Prepare for bedpostx                           #
######################################################################

# create bedpost directory inside subject and copy files to there
cd $SUBDIR
mkdir -p bedpostdata
cp "018_${SUBDIR}_b0b1000b2000.nii.gz" ./bedpostdata/data.nii.gz
cp "017_${SUBDIR}_clean_mask.nii.gz" ./bedpostdata/nodif_brain_mask.nii.gz
cp "018_${SUBDIR}_b0b1000b2000_bvals.bval" ./bedpostdata/bvals
cp "018_${SUBDIR}_b0b1000b2000_bvecs.bvec" ./bedpostdata/bvecs


# check if its all good. It should be.
# bedpostx_datacheck bedpostdata

########################### STEP 1 ###################################
#             Fiber orientation estimation                           #
######################################################################

# Creates the directory needed for the probablistic tractography
bedpostx ./bedpostdata -model 2 -n 4

cd $MAINDIR
