```
.
├── FeatureExtraction
│   ├── clean_features_from_artifacts.m
│   └── feature_extraction_full_falling_asleep.m
├── PostProcessing
│   ├── AddPostHoc_ToPrediction.m
│   ├── PredictionTipping_Summary.m
│   ├── Prediction_SDynamic.m
│   └── TippingEval_RawS_R2Summary.m
├── PredictionAnalysis
│   ├── BifurcationFittingAll
│   │   └── Calculating_R_Square_Crossings_Stats.m
│   ├── SVariableCalculation
│   │   ├── S_variable_calculation_on_clean.asv
│   │   └── S_variable_calculation_on_clean.m
│   └── UnseenNightsPrediction
│       ├── all_participants_mean_variable_script.asv
│       ├── all_participants_mean_variable_script.m
│       ├── fix_sleep_stats.asv
│       └── fix_sleep_stats.m
└── RawPSGReading
    ├── extractEEG_baseline_bedtime_so.m
    ├── extracting_full_falling_asleep_trajectory.m
    ├── find_sleep_onset.m
    └── parse_hyprograms.m
```
## RawPSGReading: How to Preprocess PSG?
**IMPORTANT: all the scripts below need directories to be adjusted in order to run!**
1. **Run `parse_hypnograms.m`.** It processes hypnogram files of participants, this will create an `AllParticipantsTable.mat`, containing the hypnogram data from each night of each of the participants row-wise
2. **Run `find_sleep_onset.m`.** It parses the `AllParticipantsTable.mat` to clean the hypnogram data and identify sleep onset in all the participants. 
3. **Run `extract_full_sleep_onset.m`.** It uses the `AllParticipantsTable.mat` as well as `{participantID}.edf` files. The output of this script are participant table files called `EpochedEEG_{participantID}.mat`, where each row represents one night, and contains the information about the night, corresponding EEG recording data for 5 minutes at the beginning of bedtime and first 5 minutes of sleep after sleep onset as well as pre-cleaned hypnogram and hypnogram indexes of the identified sleep onset 
4. **Run `extract_full_sleep_onset.m`.** This script uses the `EpochedEEG_{participantID}.mat` as well as `{participantID}.edf`. The output of this script are participant table files called `ArtifactFtEEGFallingAsleep_{participantID}.mat` where each row represents one night, and contains the information about the night, corresponding EEG recording data as well as pre-cleaned hypnogram and identified sleep onset. 
## FeatureExtraction: Extracting EEG features 
1. **Run `feature_extraction_full_falling_asleep.m`.**  It requires `ArtifactFtEEGFallingAsleep_{participantID}.mat` files to run and outputs an updated version of `ArtifactFtEEGFallingAsleep_{participantID}.mat` with feature-based representations of EEG.
2. **Run `feature_extraction_full_falling_asleep.m`.**  It requires `ArtifactFtEEGFallingAsleep_{participantID}.mat` files to run and outputs an updated version of `ArtifactFtEEGFallingAsleep_{participantID}.mat` with artifact-cleaned version of feature-based representations of EEG.

## PredictionAnalysis: How to use Feature-Based EEG representation for Sleep Distance calculation, bifurcation fitting and prediction? 
### Sleep Distance (S-Variable) Calculation 
**Run `S_variable_calculation_on_clean.m`.**  It requires `ArtifactFtEEGFallingAsleep_{participantID}.mat` files to run and outputs an updated version of `ArtifactFtEEGFallingAsleep_{participantID}.mat` with calculated true S-Variable.

### Bifurcation Fitting on true S-Variable
**Run `Calculating_R_Square_Crossing_Stats.m`.**  It requires  `ArtifactFtEEGFallingAsleep_{participantID}.mat` files with calculated true S-variables to run and outputs `FittingStatsTable.mat` which contains statistics about the fitted bifurcations (i.e. goodness of the fit, detected critical point, latency of critical points from true sleep onset)

### Sleep Onset Prediction on Unseen Nights 
1) **Run `all_participants_mean_variable_script.m`.**  It requires  `ArtifactFtEEGFallingAsleep_{participantID}.mat` files with calculated true S-variables to run and outputs `WithTestNightsNumbersSummaryTable.mat` which contains predicted S-variable for different randomly selected unseen nights and different number of nights used for modelling. Also contains a cosine similarity metric to assess how close the predicted S-variable is to the true S-variable for a given night 
2) **Run `fix_sleep_stats.m`.**  It requires `WithTestNightsNumbersSummaryTable.mat`  and outputs `uniqueTable.mat`. This file is a subset of  `WithTestNightsNumbersSummaryTable.mat` with a chosen number of nights for modelling. It contains additional columns which have data about the bifurcation fitted to the predicted S-variable, statistics about predicted critical point (predicted transition to sleep), false alarm rate as well as sleep characteristics derived from the participant's sleep staging data. 