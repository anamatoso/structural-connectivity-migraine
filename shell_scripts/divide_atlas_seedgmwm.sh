set -e
MAINDIR=$(pwd)
display_usage() {
	echo "$(basename $0) [atlas.mif] [dir_of_atlas_rois]"
	# echo "This script uses FSL to analyze diffusion data. It requires X arguments: 
	# 	1) Atlas
	}

	 if [ $# -le 1 ] # if there are 1 argument or less
	  then
	  	display_usage
	 	exit 1
	  fi

ATLAS=$1 
DIR=$2
mkdir -p "${MAINDIR}/${DIR}/atlas_rois_seedwmgm"
ATLASDIR="${MAINDIR}/${DIR}/atlas_rois_seedwmgm"

rm -f "${MAINDIR}/${DIR}/list_rois_gmwminterface.txt"

for ((i = 1 ; i <= 116 ; i++)); do
	#divide atlas into regions of interest
    mrcalc ${ATLAS} $i -eq "${ATLASDIR}/atlas_roi${i}.mif" -force # 1 if =i, 0 otherwise
	mrconvert "${ATLASDIR}/atlas_roi${i}.mif" "${ATLASDIR}/atlas_roi${i}.nii.gz" -force
	rm "${ATLASDIR}/atlas_roi${i}.mif"

	#list the files in a txt
	printf "${ATLASDIR}/atlas_roi${i}.nii.gz\n" >> "${DIR}/list_rois_gmwminterface.txt"
done
