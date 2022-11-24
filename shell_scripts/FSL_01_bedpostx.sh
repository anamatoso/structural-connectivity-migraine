#!/bin/bash

set -e #Stop on errors

display_usage() {
	echo "This script uses FSL's bedpostx to create the FODfs."
	echo "Usage: $(basename $0) [Subject and type of session]"
	echo "t requires 1 argument: the subject directory. Example: sub-control019_ses-midcycle"
	}

	if [ $# -le 0 ] # if there are 0 argument or less
	then
		display_usage
		exit 1
	fi

# Have to be in data directory
MAINDIR=$(pwd)

SUBDIR=$1 #example name: sub-control019_ses-midcycle

cd $SUBDIR

if [ ! -d "bedpostdata.bedpostX" ] # make sure bedpostx has not been run
then 
	########################### STEP 1 ###################################
	#             		   Prepare Directories                           #
	######################################################################

	# create bedpost directory inside subject and copy files there
	mkdir -p bedpostdata
	cp *_${SUBDIR}_clean.nii.gz ./bedpostdata/data.nii.gz
	cp *_${SUBDIR}_clean_mask.nii.gz ./bedpostdata/nodif_brain_mask.nii.gz
	cp "/home/amatoso/tese/data/bvals_132dir.bval" ./bedpostdata/bvals
	cp *_${SUBDIR}_eddy.eddy_rotated_bvecs.bvec ./bedpostdata/bvecs

	# check if its all good. It should be.
	# bedpostx_datacheck bedpostdata

	########################### STEP 1 ###################################
	#             		Fiber orientation estimation                     #
	######################################################################

	# Run bedpostx
	bedpostx ./bedpostdata -model 2 -n 4
	
else 
	echo "The subject ${SUBDIR} has already the bedpost directory."
fi

cd $MAINDIR