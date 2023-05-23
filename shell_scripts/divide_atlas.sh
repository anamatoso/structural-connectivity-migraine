#!/bin/bash

set -e # Stop on errors


display_usage() {
	echo ""
	tput bold 
	echo "Description"
	tput sgr0
	echo ""
	echo "This script takes the atlas file and creates a directory in which binary masks for each ROI of the atlas will be placed. It will also create a text file with the list of the files in this directory to ease probtrackx use."
	tput bold 
	echo ""
	echo "Usage"
	tput sgr0
	echo ""
	echo "./$(basename $0) [atlas.mif] [dir_of_atlas_rois]"
	echo ""
	echo "It requires 2 arguments: the atlas and the directory (that will be located inside the current directory, pwd) where the files will be placed."
	echo ""
	}

	 if [ $# -le 1 ] # if there is 1 argument or less
	  then
	  	display_usage
	 	exit 1
	  fi

MAINDIR=$(pwd)
ATLAS=$1 
DIR=$2

# Create directory where the files will be placed
ATLASDIR="${MAINDIR}/data/${DIR}/atlas_rois"
mkdir -p "${ATLASDIR}"

# Remove text if the file exists
rm -f "${MAINDIR}/data/${DIR}/list_rois.txt" 

for ((i = 1 ; i <= 116 ; i++)); do

	#divide atlas into regions of interest in nii.gz format
    mrcalc ${ATLAS} $i -eq "${ATLASDIR}/atlas_roi${i}.mif" -force # 1 if =i, 0 otherwise
	mrconvert "${ATLASDIR}/atlas_roi${i}.mif" "${ATLASDIR}/atlas_roi${i}.nii.gz" -force
	rm "${ATLASDIR}/atlas_roi${i}.mif"

	#list the files in a txt for them to be acesses in probtrackx
	printf "${ATLASDIR}/atlas_roi${i}.nii.gz\n" >> "${DIR}/list_rois.txt"
done
