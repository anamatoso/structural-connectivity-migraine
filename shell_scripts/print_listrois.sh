set -e
MAINDIR=$(pwd)
display_usage() {
	echo "$(basename $0) [atlas.mif] [dir_of_atlas_rois]"
	# echo "This script uses FSL to analyze diffusion data. It requires X arguments: 
	# 	1) Atlas
	}

	 if [ $# -le 0 ] # if there are 1 argument or less
	  then
	  	display_usage
	 	exit 1
	  fi

DIR=$1
MAINDIR2="/strombolihome/amatoso/tese/data"
rm -f "${MAINDIR}/${DIR}/list_rois.txt"

for ((i = 1 ; i <= 116 ; i++)); do
	#list the files in a txt
	printf "${MAINDIR2}/${DIR}/atlas_rois/atlas_roi${i}.nii.gz\n" >> "${DIR}/list_rois.txt"
done
