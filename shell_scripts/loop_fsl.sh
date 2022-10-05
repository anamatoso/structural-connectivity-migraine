
# MUST BE IN DATA DIRECTORY
set -e

SUBS_midcycle="sub-control019_ses-midcycle sub-control020_ses-midcycle sub-control025_ses-midcycle sub-control026_ses-midcycle sub-control027_ses-midcycle sub-control028_ses-midcycle sub-control029_ses-midcycle sub-control030_ses-midcycle sub-control031_ses-midcycle sub-control033_ses-midcycle sub-control044_ses-midcycle sub-control046_ses-midcycle sub-control048_ses-midcycle sub-control049_ses-midcycle sub-control051_ses-midcycle"
SUBS_interictal="sub-patient002_ses-interictal sub-patient003_ses-interictal sub-patient005_ses-interictal sub-patient006_ses-interictal sub-patient007_ses-interictal sub-patient008_ses-interictal sub-patient009_ses-interictal sub-patient012_ses-interictal sub-patient013_ses-interictal sub-patient034_ses-interictal sub-patient038_ses-interictal sub-patient041_ses-interictal sub-patient043_ses-interictal sub-patient045_ses-interictal"

SUBS_premenstrual="sub-control019_ses-premenstrual sub-control020_ses-premenstrual sub-control025_ses-premenstrual sub-control026_ses-premenstrual sub-control027_ses-premenstrual sub-control028_ses-premenstrual sub-control029_ses-premenstrual sub-control030_ses-premenstrual sub-control031_ses-premenstrual sub-control033_ses-premenstrual sub-control035_ses-premenstrual sub-control044_ses-premenstrual sub-control046_ses-premenstrual sub-control048_ses-premenstrual sub-control049_ses-premenstrual sub-control051_ses-premenstrual"
SUBS_ictal="sub-patient003_ses-ictal sub-patient005_ses-ictal sub-patient006_ses-ictal sub-patient013_ses-ictal sub-patient038_ses-ictal sub-patient034_ses-ictal sub-patient041_ses-ictal sub-patient045_ses-ictal sub-patient052_ses-ictal"


all="$SUBS_midcycle $SUBS_interictal $SUBS_premenstrual_resto $SUBS_ictal_resto"

for i in $all; do
    if [ -d "/home/amatoso/tese/data/${i}/bedpostdata.bedpostX" ] && [ ! -d "/strombolihome/amatoso/tese/data/${i}/bedpostdata.bedpostX" ]; then
        cp -r "/home/amatoso/tese/data/${i}/bedpostdata.bedpostX" "/strombolihome/amatoso/tese/data/${i}"
    fi
    #./FSL_01_bedpostx.sh $i
    ./FSL_02_tractography.sh $i
done


