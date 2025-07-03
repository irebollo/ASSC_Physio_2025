Code for the ASSC27 tutorial "Methods for Analyzing Brain-Body Interactions in Consciousness Research - Part 2: Analyzing Brain-Heart Interactions." by Marie Loescher.

The code needs FieldTrip and was tested on FieldTrip 20221212 and MATLAB 2022b. Other functions called by the scripts are contained in the 'Functions' folder.

If you also have the Data folder, you can run the scripts T1 to T4 on the example subject sub-26 provided in the folder, and on all subjects' GLM results. To run the code:

- Open TUTORIAL_setPaths to adapt paths to your system. This function will be called by all others to define paths.

- T1_GLM_HERs computes HERs and runs one GLM per subject, we run it only on sub-26 (provided in data folder)

- T2_GroupStats_HERs runs group statistics on the GLM results (ie beta space-time-series) of all subjects (provided in Data folder)

- T3_permutationControl_computeHERs shows you how to shuffle R-peaks for a permutation control analysis. We run it only on sub-26 (provided in data folder)

- T4_permutationControl_GroupStats shows you how to run the permutation control analysis. We can run the code up to line 174 on our example sub-26 (provided in data folder), but for the rest with group statistics we would need all subjects' data, so it doesn't run.