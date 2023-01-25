# Structural Connectivity Changes in Episodic Migraine

This repository consists of the code used in order to detect structural alterations in episodic migraine using graph theory metrics. It was done in the scope of my Master Thesis in Biomedical Engineering.

It is divided into 3 main files: `main_analysis.m`, `main_comparisons.m` and `ISMRM23.m`.

The [main_analysis.m](https://github.com/anamatoso/structural-connectivity-migraine/blob/main/matlab_scripts/main_analysis.m) file is responsible for the analysis between groups. It is dependent on several functions that were developed by me. 
They can be divided in 6 groups/objectives according to the following images:

![](https://user-images.githubusercontent.com/78906907/214705744-1ec1df5a-55ac-475a-892b-742f8e491cd2.png) <span style="font-weight:normal">Loading data</span> | ![](https://user-images.githubusercontent.com/78906907/214707275-552f4e31-be33-4a00-ad62-fc076e629f0c.png) <span style="font-weight:normal">Altering Matrices</span>
:-------------------------:|:-------------------------:
  ![](https://user-images.githubusercontent.com/78906907/214707279-007df7df-e9f5-40e8-abe1-d0e5e517a7b3.png) Calculating Connectivity Metrics | ![](https://user-images.githubusercontent.com/78906907/214707281-de1cb200-1b16-4646-a1ed-3a5d4d7fb72b.png) Visualization of Results using Boxplots
  ![](https://user-images.githubusercontent.com/78906907/214707284-82ee5342-a956-4002-9cb2-42d3b3bdc478.png) Visualization of Results using BrainNet | ![](https://user-images.githubusercontent.com/78906907/214707286-512289de-c7ff-4410-aeda-0dce3e086312.png) Statistical Analysis

Additionally, the [main_comparisons.m](https://github.com/anamatoso/structural-connectivity-migraine/blob/main/matlab_scripts/main_comparisons.m) was used to compare values between normalizations. Besides the functions mentioned before, the following were also used.

<img width="500" src="https://user-images.githubusercontent.com/78906907/214710508-390ca3e2-f284-4928-8592-e1ce46e87773.png">

Finally, the [ISMRM23.m](https://github.com/anamatoso/structural-connectivity-migraine/blob/main/matlab_scripts/ISMRM23.m) file is very similar the the main analysis however, the metrics used were calculated using a null model as a normalization. Thus, the following functions were used (or adaptedfrom others).

<img width="500" alt="Captura de ecrã 2023-01-25, às 22 58 47" src="https://user-images.githubusercontent.com/78906907/214711048-d18c6ddc-557b-40c8-be56-0666bb240cf9.png">

Furthermore, if you need any further explanation, do not hesitate to contact me!
