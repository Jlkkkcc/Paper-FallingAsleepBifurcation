# Paper-FallingAsleepBifurcation
 The scripts and final data repository for paper: "Falling Asleep displays a predictable bifurcation dynamic".
 
 **Potential competing interests**:
JL and NG have a patent application filed (Patent Application No. 2509587.8) for the methods and computational framework here.

**Paper DOI**:

## Contents

- [Sleep datasets used](#sleep-datasets-used)
- [Requirements check](#requirements-check)
- [Basic introduction of this repository](#basic-introduction-of-this-repository)
- [Revision Notes May-2025](#revision-notes-may-2025)
- [Folder details](#folder-details)


## Sleep datasets used
**The raw sleep datasets**:

1. Multi-Ethnic Study of Atherosclerosis (MESA); Open-source dataset that can be obtained from National Sleep Research Resources (NSRR); https://sleepdata.org/datasets/mesa 

2. The cohort 2 dataset was collected by the authors, and it can be made available under sufficient material transfer agreements (MTA). Please contact the authors.


## Requirements check
The code was written and is working in MATLAB R2023a version and Python 3.8, with the additional functions/packages all working under this environment and tested in the Ubuntu system. If the code does not work or run on your side, check firstly your software versions (including external software not included in the repository, please download and install them on your side); Contact the authors if you have done every essential check. Please also ensure that you have installed all relevant MATLAB internal toolboxes (if any are missing).

Sometimes the external software has been updated, so necessary changes are necessary to make the code run. We have tried to mark them as far as we knew; however, it might not be comprehensive. 

### External Packages required

1. **CATCH-22**; Please download and install following their instructions to compile the C codes in MATLAB [here](https://github.com/DynamicsAndNeuralSystems/catch22).

   **Notes:** This package has been updated from the version we used (27 features) to 29 features. Excluding the last two features, or modify the codes if you want to include the last two features. The version used in the manuscript analysis is v0.1.0. Please also note, _this toolbox update would result in quantitative differences in the feature outputs, as well as different sequences of feature outputs (compared to the manuscript results)_, but should not affect bifurcation dynamic findings.
   
   Typical installation time: 1-2 minutes
2. Python package for Functional principal component analysis (FPCA) [here](https://fda.readthedocs.io/en/latest/).

	The package dependencies for the FPCA analysis that generated the results in the manuscript: Python 3.8, Scipy (Version 1.6.2), Scikit-fda (Version 0.7)

	**Notes:** This package has been updated if you use Python 3.10, and the code given might return errors running Regularised FPCA analysis.
	
	Typical installation time: 1-minute
3. **EntRate** -- Entropy rate estimators for neuroscience (normalised Lempel-Ziv complexity), download and install [here](https://github.com/pmediano/EntRate). Similarly, please follow the instructions to compile the C code.
	
	Typical installation time: <1-minute
4. **FOOOF**: a python toolbox for evaluating the aperiodic (1/f) component of the EEG power spectrum; Please follow instructions to download and install [here](https://fooof-tools.github.io/fooof/#).
	
	**Notes:** Despite that this is a Python toolbox, this was used in MATLAB, with a Python interface setup. Please ensure that you have the corresponding Python interface set up and functioning before attempting to run the code successfully.
	
	The Python interface used in the paper: Python 3.10 with basic packages installed (see requirements on the link).
	
	Typical installation time: 1-minute


## Basic introduction of this repository

Before you start to run analyses or plotting, _make sure that you add the entire folder to the path in MATLAB_. The repository contains four folders:

1. **"*Core Algorithms and Examples*"**: This folder contains the two core algorithms in this paper: feature and sleep distance evaluation, the bifurcation function fitting, and the prediction of the tipping point. The folder also includes some example data for exploration on your own.

2. **"*Process Codes (NoData)*"**: This folder contains the process codes to read, process, and analyse the raw sleep data (_everything from scratch_). However, the raw sleep data (PSG) should be obtained as specified previously; no raw data has been deposited here.

3. **"*Final Plot codes*"**: This folder includes all the codes for generating the main figures as well as Supplementary figures (excluding demographics) of the manuscript.

4. **"*All Final Data for plots*"**: This folder only contains all the data necessary for generating plots (Used for codes in the 'Final plots codes' folder). No raw data has been included. 

5. **"*Revision analysis (add features)*"**: This folder contains the new codes for the extraction and patching (to the original data structure) of three new EEG features: Sigma band power, Lempel-Ziv complexity, and the aperiodic (1/f) component.

#### **Instructions for exploration:**
If you would like to see how the core algorithms in the paper work, please go to the folder **"*Core Algorithms and Examples*"**. If you would like to run through analysis from scratch (from raw PSG data to the end results), you can go to **"*Process Codes (NoData)*"**. If you would like to see figures replicated, please go to **"*Final Plot codes*"**. The other folders document the data and revision analysis codes. 

## Revision Notes May-2025

There were two types of analysis during revision. The biggest revision was adding features, where the added feature codes are inside a specific folder, "Revision analysis (add features)". The other analyses were either revising original results (replacing previous files) or were new results directly added to the analysis folders.

## Folder details

#### Core Algorithms and Examples

The folder provides three examples of core algorithms, organised in separate folders. Before running the code, remember to install the relevant toolboxes and add the entire repository to the directory.

1. Feature extraction and sleep distance calculation in folder '_Feature and Sleep distance computation_'.

	The "`Feature_SDist_Eval.m`" provides an exemplary code to extract EEG features from EEG epochs (6 s) for the falling asleep process (cut from 30 minutes before sleep onset to 10 minutes after), and the evaluation of the sleep distance metric. This same sleep distance time series is used for the bifurcation fitting example (see next).
	
	***Important notice:*** This code is not a replication of the full computational pipeline but only the key computational aspects above. The other processing steps (e.g., artefact removal) are not included in this example; the EEG data provided has already been processed, and we here only included one EEG channel to shorten the runtime. Please do follow standard guidelines and your own demands to pre-process your EEG data before the feature extraction step.

	The typical runtime for this code will be *~10 minutes* on a "normal" computer. Notice that the feature extraction is very time-consuming, therefore, we also included the pre-evaluated features for evaluation of sleep distance (*<1 minute*).

2. Bifurcation fitting examples in folder '_Examples Bif Fitting_'. 

	The "`BifurcationFitting_Example.m`" gives an example of data for you to explore how different initial parameters could lead to different fits; It then optimises the fit based on our in-house developed functions. After fitting, the tipping point was evaluated. The key algorithm functions and required data are included inside.

   The typical runtime for this code will be *<1 minute* on a "normal" computer for bifurcation fitting, and *1-2 minutes* for tipping point evaluation.

3. Tipping point prediction example in folder '_Example Prediction Tipping point_'. 

	The "`Plot_PredExample.m`" shows one complete participant, with one training night and seven testing nights. It shows and plots the prediction of tipping points compared to the post-hoc model results. The key algorithm functions (in-house developed) and required data are included inside.

	The typical runtime for this code will be *1-3 minutes* on a "normal" computer.


#### Final Plot codes & All Final Data for plots

These two folders are well organised in the order of the figure. For Figure 3, the plot codes are separated given the complexity of different features. Make sure that you added all functions and data to the path.

In 'Figure3' Folder:

1. "`FPCA_Features.py`": The FPCA analysis and the plots in one file (Figure 3a and 3d).
2. "`Plot_IndividualFts.m`": The plots of individual feature time series dynamics (Figure 3b and 3e).
3. "`SleepDistance\_Features_Plot.m`": The plots where features are plotted against the sleep distance dynamics (Figure 3c and 3f).

#### Process Codes (From Raw PSG data)

The process codes are divided further into 3 subfolders, containing all the functions used ("*All functions*"), the analysis in the MESA dataset ("*MESA analysis*") and the analysis in cohort 2 ("*Cohort2 Analysis*"). 

There are a few externally developed packages in MATLAB (which do not require installation) in the "*All functions*" folder:

* [ShadedErrorBar](https://uk.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar)
* [ViolinPlot](https://github.com/bastibe/Violinplot-Matlab)
* [DistinguishableColors](https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors#functions_tab)
* [For PSG (EDF format) Reading](https://www.edfplus.info/downloads/index.html)

The other subfolders within contain our self-developed in-house algorithms and codes for the relevant analysis required.

> The codes are commented within. If you have any questions about using or about detailed variable explanations, please get in touch with us.

----
#### Analysis codes of the MESA dataset

The analysis scripts are introduced in a sequential order below:

* ***RawPSG_Loading***: This folder contains the MATLAB codes for reading the raw PSG files (*`Data_load.m`*), as well as participant exclusion checks (*`Demo_Sleep_Checks.m`*), and demographic analysis (*`Demo_Checks.m`*). This step will return a list of Matlab data files for feature extraction in the second step.

* ***Feature Extractions***: This folder contains the code for feature extraction (*`fteval_allsbj.m`*), using the data coming from the first step.

* ***FeatureSpaceAnalysis***: This folder contains the main analysis code to generate results for Figure 2-3 in the manuscript.
	1. *Preparation*: This folder contains the script (*`FeatureDynam_Prep.m`*) to prepare the features from the second step towards further group-level feature-space analysis. They include artefact checks on EEG epochs, as well as global z-score normalisation.
	2. *MainVariables\_Group_Figure2*: Group-level analysis of the feature space variables including:
		* Sleep distance `s(t)`: *_`SleepDistanceS.m`_*
		* State velocity `v(t)`: *_`StateVelocityVt.m`_*
		* Sleep velocity `v_s(t)`: *_`SleepVelocityVs.m`_*
		* Bifurcation function fitting of the sleep distance: *_`BifFitting_SleepDistance.m`_*
	
		These codes generate the results for parts of Figures 2 and 3.
	
	3. *FeatureAnalyis_Figure3*: Group-level analysis of individual feature dynamics:
	   * FPCA analysis on the group: *_`FPCA_Features_Grp.py`_* and *_`FPCAscores_anal.m`_*
	   * FPCA analysis conducted on the individual level: *_`TSPrep_IndivFPCA.m`_* and *_`Individual_FPCA.py`_*
	   
	   These codes generate the results for part of Figure 3. 
	   
	4. *IndividualAnalysis_Figure2*: Analysis of the individual-level sleep-distance dynamics (for results in Figure 2f-i):
		* *_`TSPrep_EachIndividual.m`_*: this code generates the Sleep distance s(t) dynamics for each individual (With some further exclusions due to artefact checks);
		* Early warning signal for Critical slowing down analysis: *_`EWS_Evaluations.m`_* and *_`EWSAnalysis_Stats.m`_*.
		* Bifurcation function fitting on individual-level dynamics analysis: *_`Individual_BifFitting.m`_* and *_`Individual_BifFitting_Analysis.m`_*. The codes include all the statistical analysis and correlation analysis results described in the manuscript.

	5. *BifAnalysis\_IndividualFeatures\_Rev*: This folder contains analysis codes for revision supplements, where the bifurcation dynamics were tested on individual features (the nine FPC1 representative features). The code *_`Individual_BifFitting_perftLoop.m`_* prepares the feature time series for bifurcation fitting, and *_`TSPrep_IndivFt_IndivSbj_Loopfts.m`_* runs the bifurcation fitting and saves the results.

To see the figure results, go to the figures folder and visualise them. Note that the many trivial analysis codes (such as the summary of certain statistics like mean or standard deviations, or quick test of correlations, as reported in the manuscript) are not all included; the codes here are aiming to cover the main results and key methodologies used.


-----
#### Analysis for Cohort 2

For the process codes (computation from the raw data), please refer to the documentation in "*Cohort2 Analysis*" folder: "`Cohort2 Code Description.md`"; Only the post-processing codes for generating figures are described below.

***PostProcessing***: This subfolder post-processes all results generated previously for displaying in figures, as well as a statistics summary, using the data from the process codes.

* 	Computing the tipping points from the post-hoc bifurcation function fits for the Sleep distance dynamics from all nights (*_`TippingEval_RawS_R2Summary.m`_*)
* 	Adding post-hoc bifurcation function fit tipping point times toward the prediction results (*_`AddPostHoc_ToPrediction.m`_*)
* 	Summary analysis of the prediction results of the sleep distance dynamics (compared to ground truth); The accuracy was evaluated using the cosine similarity scores; (*_`Prediction_SDynamic.m`_*)
* 	Summary analysis of the tipping point prediction results (*_`PredictionTipping_Summary.m`_*)

To see the analysis results, go to the figures folder and visualise them.

	
















