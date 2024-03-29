#!/bin/bash

set -e #Stop on errors

display_usage() {
	tput bold 
	echo ""
	echo "Description"
	tput sgr0
	echo ""
	echo "This script uses FSL's bedpostx to create the FODfs. It will create the bedpostdata.bedpostX directory inside the input directory."
	tput bold 
	echo ""
	echo "Usage"
	tput sgr0
	echo ""
	echo "./$(basename $0) [Subject and type of session]"
	echo ""
	echo "It requires 1 argument: the subject directory. Example: sub-control019_ses-midcycle"
	}

	if [ $# -le 0 ] # if there are no arguments
	then
		display_usage
		exit 1
	fi

# Have to be in data directory
MAINDIR=$(pwd)

SUBDIR=$1
cd "data/${SUBDIR}"

if [ ! -d "bedpostdata.bedpostX" ] # make sure bedpostx has not been run
then 
	########################### STEP 1 ###################################
	#             		   Prepare Directories                           #
	######################################################################

	# create bedpost directory inside subject and copy files there (with the names required by bedpostx)
	# this directory should have the data (dwi), the brain mask, the bvals and bvecs
	mkdir -p bedpostdata
	cp *_${SUBDIR}_clean.nii.gz ./bedpostdata/data.nii.gz
	cp *_${SUBDIR}_clean_mask.nii.gz ./bedpostdata/nodif_brain_mask.nii.gz
	cp "${MAINDIR}/bvals_132dir.bval" ./bedpostdata/bvals
	cp *_${SUBDIR}_eddy.eddy_rotated_bvecs.bvec ./bedpostdata/bvecs

	# To check if its all good use the following command: bedpostx_datacheck bedpostdata

	########################### STEP 2 ###################################
	#             		Fiber orientation estimation                     #
	######################################################################

	# Run bedpostx
	bedpostx ./bedpostdata -model 2 -n 4 # maximum of 4 fibers per voxel
	
else 
	echo "The subject ${SUBDIR} has already gone through bedpostx."
fi

cd $MAINDIR