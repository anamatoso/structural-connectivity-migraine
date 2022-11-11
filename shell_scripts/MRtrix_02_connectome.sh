#!/bin/bash

MAINDIR=$(pwd)
display_usage() {
	echo "$(basename $0) [Subject] [type of session]"
	# echo "This script uses FSL to analyze diffusion data. It requires X arguments: 
	# 	1) Subject DWI Directory (string)- has to have data.nii.gz (...b0b1b2000.nii.gz), nodif_brain_mask.nii.gz, bvecs (...b0b1000b2000), bvals
	#	example: "sub-control019"
	#	2) Type of session. Example: ses-midcycle
	}

	if [ $# -le 1 ] # if there are 1 argument or less
	then
		display_usage
		exit 1
	fi


SUB=$1 #example name: sub-control019
SES=$2

########################### STEP 1 ###################################
#             		  Get data from directories                      #
######################################################################

SUBDIR="${MAINDIR}/${SUB}_${SES}" #example name: sub-control019_ses-midcycle

ANATDIR="${MAINDIR}/${SUB}"
ANAT="${ANATDIR}/${SUB}_restored-MPRAGE_brain.nii.gz"

DWI="${MAINDIR}/${SUB}_${SES}/*b0b1000b2000.nii.gz"
DWIMASK="${MAINDIR}/${SUB}_${SES}/*clean_mask.nii.gz"
BVEC="${MAINDIR}/${SUB}_${SES}/*b0b1000b2000_bvecs.bvec"
BVAL="${MAINDIR}/${SUB}_${SES}/*b0b1000b2000_bvals.bval"

