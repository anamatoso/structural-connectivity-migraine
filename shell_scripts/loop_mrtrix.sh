
# MUST BE IN DATA DIRECTORY
set -e

SUBS_midcycleOG="sub-control019_ses-midcycle sub-control020_ses-midcycle sub-control025_ses-midcycle sub-control026_ses-midcycle sub-control027_ses-midcycle sub-control028_ses-midcycle sub-control029_ses-midcycle sub-control030_ses-midcycle sub-control031_ses-midcycle sub-control033_ses-midcycle sub-control044_ses-midcycle sub-control046_ses-midcycle sub-control048_ses-midcycle sub-control049_ses-midcycle"
SUBS_interictal="sub-patient002_ses-interictal sub-patient003_ses-interictal sub-patient005_ses-interictal sub-patient006_ses-interictal sub-patient007_ses-interictal sub-patient008_ses-interictal sub-patient009_ses-interictal sub-patient012_ses-interictal sub-patient013_ses-interictal sub-patient034_ses-interictal sub-patient038_ses-interictal sub-patient041_ses-interictal sub-patient043_ses-interictal sub-patient045_ses-interictal"
SUBS_midcycle="sub-control025_ses-midcycle sub-control026_ses-midcycle sub-control027_ses-midcycle sub-control028_ses-midcycle sub-control029_ses-midcycle sub-control030_ses-midcycle sub-control031_ses-midcycle sub-control033_ses-midcycle sub-control044_ses-midcycle sub-control046_ses-midcycle sub-control048_ses-midcycle sub-control049_ses-midcycle"


for SUBDIR in $SUBS_midcycle; do
	./MRTrix_Script_bvals2_newatlas.sh $SUBDIR
done

for SUBDIR in $SUBS_interictal; do
	./MRTrix_Script_bvals2_newatlas.sh $SUBDIR
done


