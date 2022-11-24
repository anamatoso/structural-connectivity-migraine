#!/bin/bash

set -e #Stop on errors


display_usage() {
	echo "This script takes the atlas files and creates files which are binary masks that correspond to each ROI of the atlas."
	echo "Usage: $(basename $0) [atlas.mif] [dir_of_atlas_rois]"
	echo "It requires 2 arguments: the atlas and the directory where the files will be places. Example: sub-control019_ses-midcycle"
	}

	 if [ $# -le 1 ] # if there are 1 argument or less
	  then
	  	display_usage
	 	exit 1
	  fi

MAINDIR=$(pwd)
ATLAS=$1 
DIR=$2

# Create directory where the files will be placed
mkdir -p "${MAINDIR}/${DIR}/atlas_rois"
ATLASDIR="${MAINDIR}/${DIR}/atlas_rois"

rm -f "${MAINDIR}/${DIR}/list_rois.txt" # Replace if the file exists

for ((i = 1 ; i <= 116 ; i++)); do
	#divide atlas into regions of interest
    	mrcalc ${ATLAS} $i -eq "${ATLASDIR}/atlas_roi${i}.mif" -force # 1 if =i, 0 otherwise
	mrconvert "${ATLASDIR}/atlas_roi${i}.mif" "${ATLASDIR}/atlas_roi${i}.nii.gz" -force
	rm "${ATLASDIR}/atlas_roi${i}.mif"

	#list the files in a txt
	printf "${ATLASDIR}/atlas_roi${i}.nii.gz\n" >> "${DIR}/list_rois.txt"
done
