
# MUST BE IN DATA DIRECTORY

# Controls have midcycle and premenstrual
SUB="sub-control019 sub-control020"
SES="ses-midcycle ses-premenstrual"


# for i in $SUB; do
#     for j in $SES; do
#         #MRTrix
#         #./MRtrix_Script.sh $i $j
#         #./FSL_01_bedpostx.sh ""
#         #FSL bedpostx
#     done
# done

SUBS="sub-control019_ses-premenstrual sub-control020_ses-midcycle sub-control020_ses-premenstrual sub-control025_ses-midcycle sub-control025_ses-premenstrual sub-control026_ses-midcycle sub-control026_ses-premenstrual sub-control027_ses-midcycle sub-control027_ses-premenstrual sub-control028_ses-premenstrual sub-control029_ses-midcycle sub-control029_ses-premenstrual sub-control030_ses-midcycle sub-control030_ses-premenstrual sub-control031_ses-midcycle sub-control031_ses-premenstrual sub-control033_ses-midcycle sub-control033_ses-premenstrual sub-control035_ses-premenstrual sub-patient002_ses-interictal sub-patient003_ses-ictal sub-patient003_ses-interictal sub-patient005_ses-ictal sub-patient005_ses-interictal sub-patient006_ses-ictal sub-patient006_ses-interictal sub-patient007_ses-interictal sub-patient008_ses-interictal sub-patient009_ses-interictal sub-patient012_ses-interictal sub-patient013_ses-ictal sub-patient034_ses-interictal sub-patient038_ses-ictal sub-patient038_ses-interictal"
SUBS_midcycle="sub-control026_ses-midcycle sub-control027_ses-midcycle sub-control029_ses-midcycle sub-control030_ses-midcycle sub-control031_ses-midcycle sub-control033_ses-midcycle"

for i in $SUBS_midcycle; do
    ./FSL_01_bedpostx.sh $i
done