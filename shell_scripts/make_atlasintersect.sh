
# MUST BE IN DATA DIRECTORY
set -e
MAINDIR=$(pwd)
mrconvert AAL116_MNI152-2mm.nii.gz AAL116_MNI152-2mm.mif -force -quiet
cp AAL116_MNI152-2mm.mif AAL116_intersect.mif

SUBS_midcycle="sub-control019_ses-midcycle sub-control020_ses-midcycle sub-control025_ses-midcycle sub-control026_ses-midcycle sub-control027_ses-midcycle sub-control029_ses-midcycle sub-control030_ses-midcycle sub-control031_ses-midcycle sub-control033_ses-midcycle sub-control044_ses-midcycle sub-control046_ses-midcycle sub-control049_ses-midcycle"
SUBS_interictal="sub-patient002_ses-interictal sub-patient003_ses-interictal sub-patient005_ses-interictal sub-patient006_ses-interictal sub-patient007_ses-interictal sub-patient008_ses-interictal sub-patient009_ses-interictal sub-patient012_ses-interictal sub-patient013_ses-interictal sub-patient034_ses-interictal sub-patient038_ses-interictal sub-patient041_ses-interictal sub-patient043_ses-interictal sub-patient045_ses-interictal"
SUBS_premenstrual="sub-control019_ses-premenstrual sub-control020_ses-premenstrual sub-control025_ses-premenstrual sub-control026_ses-premenstrual sub-control027_ses-premenstrual sub-control028_ses-premenstrual sub-control029_ses-premenstrual sub-control030_ses-premenstrual sub-control031_ses-premenstrual sub-control033_ses-premenstrual sub-control035_ses-premenstrual sub-control046_ses-premenstrual sub-control048_ses-premenstrual"
SUBS_ictal="sub-patient003_ses-ictal sub-patient005_ses-ictal sub-patient006_ses-ictal sub-patient013_ses-ictal sub-patient038_ses-ictal sub-patient034_ses-ictal sub-patient041_ses-ictal"

alls_SUBS="sub-control019_ses-midcycle sub-control020_ses-midcycle sub-control025_ses-midcycle sub-control026_ses-midcycle sub-control027_ses-midcycle sub-control029_ses-midcycle sub-control030_ses-midcycle sub-control031_ses-midcycle sub-control033_ses-midcycle sub-control044_ses-midcycle sub-control046_ses-midcycle sub-control049_ses-midcycle sub-patient002_ses-interictal sub-patient003_ses-interictal sub-patient005_ses-interictal sub-patient006_ses-interictal sub-patient007_ses-interictal sub-patient008_ses-interictal sub-patient009_ses-interictal sub-patient012_ses-interictal sub-patient013_ses-interictal sub-patient034_ses-interictal sub-patient038_ses-interictal sub-patient041_ses-interictal sub-patient043_ses-interictal sub-patient045_ses-interictal sub-control019_ses-premenstrual sub-control020_ses-premenstrual sub-control025_ses-premenstrual sub-control026_ses-premenstrual sub-control027_ses-premenstrual sub-control028_ses-premenstrual sub-control029_ses-premenstrual sub-control030_ses-premenstrual sub-control031_ses-premenstrual sub-control033_ses-premenstrual sub-control035_ses-premenstrual sub-control046_ses-premenstrual sub-control048_ses-premenstrual sub-patient003_ses-ictal sub-patient005_ses-ictal sub-patient006_ses-ictal sub-patient013_ses-ictal sub-patient038_ses-ictal sub-patient034_ses-ictal sub-patient041_ses-ictal"
ATLAS="/home/amatoso/tese/data/AAL116_MNI152-2mm.nii.gz"

for SUBDIR in $alls_SUBS; do
    ANATDIR="${MAINDIR}/${SUBDIR:0:14}"
    cd $SUBDIR
    DWIMASK="*clean_mask.nii.gz"
    DWI="015_${SUBDIR}_clean.nii.gz"
    # # tranform from diff space to struct space
    # mrconvert $DWIMASK dwimask.mif -force -quiet
    # mrtransform dwimask.mif -linear "mrtrix_outputs_bvals2/diff2struct_mrtrix.txt" dwimask2struct.mif -force -quiet
    # mrconvert dwimask2struct.mif dwimask2struct.nii.gz -force

    # # tranform from struct space to standard space
    # applywarp -i dwimask2struct.nii.gz -r $ATLAS --out=mask_2standard --rel --warp="${ANATDIR}/reg_nonlinear_warp_T1tostandard_2mm.nii.gz"
    # mrconvert mask_2standard.nii.gz mask_2standard.mif -force -quiet

    fnirt --ref=$ATLAS --in=$DWI --cout=diff2standard.nii.gz
    applywarp -i $DWIMASK -r $ATLAS --out=mask_2standard --warp=diff2standard.nii.gz
    mrconvert mask_2standard.nii.gz mask_2standard.mif -quiet

    #rm -f mask_2standard.nii.gz
    rm -f dwimask.mif
    rm -f dwimask2struct.mif
    rm -f flirt_affine_guess.mat
    rm -f diff2standard.nii.gz
    
    cd ..
    mrcalc AAL116_intersect.mif "${MAINDIR}/${SUBDIR}/mask_2standard.mif" -mult temp.mif -quiet
    cp temp.mif AAL116_intersect.mif
    rm temp.mif
    rm -f "${MAINDIR}/${SUBDIR}/mask_2standard.mif"
    echo $SUBDIR
done

mrconvert AAL116_intersect.mif AAL116_intersect.nii.gz