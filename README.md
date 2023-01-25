# Structural Connectivity Changes in Episodic Migraine
<img width="1271" alt="banner" src="https://user-images.githubusercontent.com/78906907/214711685-eff796e6-f530-40c3-b40d-eb39f177d952.png">

This repository consists of the original code created in order to detect structural alterations in episodic migraine using graph theory metrics. It was done in the scope of my Master's Thesis in Biomedical Engineering at Instituto Superior Técnico. In this work, a comparison was made not only between controls and patients but also longitudinally between different stages of the menstrual cycle. In addition, it was also studied the impact of different software packages on the calculation of the connectivity matrix and of different normalizations of the connectivity matrix.


Firstly, from the DWI images, connectivity matrices had to be created. To do that, the bash files in shell_scripts were used. Two software packages were used and compared: MRtrix and FSL. To create a connectivity matrix using MRtrix one only needs to run the `MRTrix_Script.sh`. On the other hand, if you want to use FSL, you only need to run `FSL_01_bedpostx.sh` followed by either `FSL_02_tractography.sh` or `FSL_02_tractography_mat3.sh`. Note that to run the tractography in FSL, you first need to run `divide_atlas.sh`.

Then, with the connectivity matrices of all subjects, an analysis in MATLAB can be done using the files in matlab_scripts.
This folder has 3 main files: `main_analysis.m`, `main_comparisons.m` and `ISMRM23.m`. These files are dependent on several functions that had to be created by me whose dependencies are explained below. 

The [main_analysis.m](https://github.com/anamatoso/structural-connectivity-migraine/blob/main/matlab_scripts/main_analysis.m) file is responsible for the analysis between groups.
It is dependent on many functions that can be divided into 6 groups/objectives according to the following images:

![](https://user-images.githubusercontent.com/78906907/214705744-1ec1df5a-55ac-475a-892b-742f8e491cd2.png) <span style="font-weight:normal">Loading data</span> | ![](https://user-images.githubusercontent.com/78906907/214707275-552f4e31-be33-4a00-ad62-fc076e629f0c.png) <span style="font-weight:normal">Altering Matrices</span>
:-------------------------:|:-------------------------:
  ![](https://user-images.githubusercontent.com/78906907/214707279-007df7df-e9f5-40e8-abe1-d0e5e517a7b3.png) Calculating Connectivity Metrics | ![](https://user-images.githubusercontent.com/78906907/214707281-de1cb200-1b16-4646-a1ed-3a5d4d7fb72b.png) Visualization of Results using Boxplots
  ![](https://user-images.githubusercontent.com/78906907/214707284-82ee5342-a956-4002-9cb2-42d3b3bdc478.png) Visualization of Results using BrainNet | ![](https://user-images.githubusercontent.com/78906907/214707286-512289de-c7ff-4410-aeda-0dce3e086312.png) Statistical Analysis

Additionally, the [main_comparisons.m](https://github.com/anamatoso/structural-connectivity-migraine/blob/main/matlab_scripts/main_comparisons.m) was used to compare values between normalizations. Besides the functions mentioned before, the following were also used:

<img width="400" src="https://user-images.githubusercontent.com/78906907/214710508-390ca3e2-f284-4928-8592-e1ce46e87773.png">

Finally, the [ISMRM23.m](https://github.com/anamatoso/structural-connectivity-migraine/blob/main/matlab_scripts/ISMRM23.m) file is very similar to the main analysis however, the metrics used were calculated using a null model as a normalization. Thus, the following functions were used:

<img width="400" alt="Captura de ecrã 2023-01-25, às 22 58 47" src="https://user-images.githubusercontent.com/78906907/214711048-d18c6ddc-557b-40c8-be56-0666bb240cf9.png">

Furthermore, if you need any further explanation, do not hesitate to contact me!
