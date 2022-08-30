set -e

all_SUBSdwi="sub-control051_ses-premenstrual"
#link dwi images
for i in $all_SUBSdwi; do
 sub="/strombolihome/mig_n2treatdata/derivatives/dwi-preproc/${i}"
 BASEDIR=$(basename "$sub")
 mkdir -p $BASEDIR
 cd $BASEDIR
 for file in /strombolihome/mig_n2treatdata/derivatives/dwi-preproc/$BASEDIR/*; do
    if [ -f "$file" ] 
    then
    BASEFILE=$(basename "$file")
    ln -s $file $BASEFILE -f
    fi
 done
 cd ..
done


all_SUBS="sub-control051"

#link anatomical images
for i in $all_SUBS; do
 sub="/strombolihome/mig_n2treatdata/derivatives/anat-preproc/fsl_fnirt_reg2standard/${i}"
 BASEDIR=$(basename "$sub")
 mkdir -p $BASEDIR
 cd $BASEDIR
 for file in /strombolihome/mig_n2treatdata/derivatives/anat-preproc/fsl_fnirt_reg2standard/$BASEDIR/*; do
    if [ -f "$file" ] 
    then
    BASEFILE=$(basename "$file")
    ln -s $file $BASEFILE -f
    fi
 done
 cd ..
done

for i in $all_SUBS; do
 sub="/strombolihome/mig_n2treatdata/derivatives/anat-preproc/ants_bet/${i}" 
 BASEDIR=$(basename "$sub")
 mkdir -p $BASEDIR
 cd $BASEDIR
 for file in /strombolihome/mig_n2treatdata/derivatives/anat-preproc/ants_bet/$BASEDIR/*; do
    if [ -f "$file" ] 
    then
    BASEFILE=$(basename "$file")
    ln -s $file $BASEFILE -f
    fi 
done
 cd ..
done
