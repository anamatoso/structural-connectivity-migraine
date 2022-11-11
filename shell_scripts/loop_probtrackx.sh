
# MUST BE IN DATA DIRECTORY
set -e

SUBS="sub-control019_ses-premenstrual sub-control020_ses-midcycle sub-control020_ses-premenstrual sub-control025_ses-midcycle sub-control025_ses-premenstrual sub-control026_ses-midcycle sub-control026_ses-premenstrual sub-control027_ses-midcycle sub-control027_ses-premenstrual sub-control028_ses-premenstrual sub-control029_ses-midcycle sub-control029_ses-premenstrual sub-control030_ses-midcycle sub-control030_ses-premenstrual sub-control031_ses-midcycle sub-control031_ses-premenstrual sub-control033_ses-midcycle sub-control033_ses-premenstrual sub-control035_ses-premenstrual sub-patient002_ses-interictal sub-patient003_ses-ictal sub-patient003_ses-interictal sub-patient005_ses-ictal sub-patient005_ses-interictal sub-patient006_ses-ictal sub-patient006_ses-interictal sub-patient007_ses-interictal sub-patient008_ses-interictal sub-patient009_ses-interictal sub-patient012_ses-interictal sub-patient013_ses-ictal sub-patient034_ses-interictal sub-patient038_ses-ictal sub-patient038_ses-interictal"
SUBS_midcycle="sub-control019_ses-midcycle sub-control020_ses-midcycle sub-control025_ses-midcycle sub-control026_ses-midcycle sub-control027_ses-midcycle sub-control029_ses-midcycle sub-control030_ses-midcycle sub-control031_ses-midcycle sub-control033_ses-midcycle"
SUBS_interictal="sub-patient002_ses-interictal sub-patient003_ses-interictal sub-patient005_ses-interictal sub-patient006_ses-interictal sub-patient007_ses-interictal sub-patient008_ses-interictal sub-patient009_ses-interictal sub-patient012_ses-interictal sub-patient013_ses-interictal sub-patient034_ses-interictal sub-patient038_ses-interictal"

SUBS_midcycle_resto="sub-control044_ses-midcycle sub-control046_ses-midcycle sub-control049_ses-midcycle"



for i in $SUBS_interictal_resto; do
    ./FSL_02_tractography.sh $i
done
