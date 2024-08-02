# Paper-FallingAsleepBifurcation
 The code and final data repository for paper: "Falling Asleep is a predictable catastrophic bifurcation dynamic".

**Paper DOI**:

**The raw sleep datasets**:

1. Multi-Ethnic Study of Atherosclerosis (MESA); Open-source dataset that can be obtained from National Sleep Research Resources (NSRR); https://sleepdata.org/datasets/mesa 

2. The cohort 2 dataset was collected by authors, and it can be made available under sufficient material transfer agreements (MTA). Please contact the authors.

**Competing interests**:
The methodology and algorithm of this manuscript are in the process of being patented. If you use it make sure it was cited.

## Some notices
The code was written and working in MATLAB R2023a version and Python 3.7, with the additional functions/packages all working under this environment and tested in the Ubuntu system. If the code does not work or run on your side, check firstly your software versions (including external software not included here); Contact the authors if you have done every essential check.

Sometimes the external software has been updated so that necessary changes are necessary to make the code running. I have tried to mark them as far as I know, however might not be comprehensive. Contact us if so.

### External Packages required

1. CATCH-22; Please download and install following their instructions [here](https://github.com/DynamicsAndNeuralSystems/catch22).

   **Notes:** This package has been updated from the version we used (27 features) now to 29 features. Excluding the last two features or modify the codes if you want to include the last two features
2. Python package for Functional principal component analysis (FPCA) [here](https://fda.readthedocs.io/en/latest/).

	**Notes:** This package has been updated if you use Python 3.10 and the code given might return errors running Regularised FPCA analysis.

## Basic introduction of this repository

Before you start to run analysis or plotting, make sure that you add the entire folder to the path in MATLAB. The repository contains four folders:

1. "*Core Algorithms and Examples*": This folder contains the two core algorithms in this paper: the bifurcation function fitting, and the prediction of the tipping point. The folder also includes some example data for exploration on your own.

2. "*Process Codes (NoData)*": This folder contains the process codes to read, process, and analysis of the raw sleep data. The raw sleep data should be obtained as specified previously, no data has been given here.

3.  "*Final Plot codes*": This folder includes all the codes for generating the main figures of the manuscript.

4. "*All Final Data for plots*": This folder contains all the data necessary for generating plots only. No raw data has been included. 

## Folder details

#### Core Algorithms and Examples

Key algorithmic functions are included in the two folders inside: "Prediction Tipping Point" and "Bifurcation Function Fitting".

The two MATLAB codes show two examples (along with two ".mat" files containing the data) in using the functions.

1. The "Plot_PredExample.m" shows one complete participant, with one training night and seven testing nights. It shows and plots the prediction of tipping points compared to the post-hoc model results;

2. The "BifurcationFitting_Example.m" gives an example of data for you to explore how different initial parameters could lead to different fits; It then contains the optimisation we developed.

#### Final Plot codes & All Final Data for plots

These two folders are well organised in the order of the figure. For Figure 3, the plot codes are separated up given the complexity of different features. Make sure that you added all functions and data to path.

In Figure3 Folder:

1. "FPCA_Features.py": The FPCA analysis and the plots in one file (Figure 3A and 3D).
2. "Plot_IndividualFts.m": The plots of individual feature time series dynamics (Figure 3B and 3E).
3. "SleepDistance\_Features_Plot.m": The plots where features are plotted against the sleep distance dynamics (Figure 3C and 3F).
4. "GC_Plot.m": the Granger-causality test results (Figure 3G).

#### Process Codes (NoData)

The process codes are divided further into 3 subfolders, containing all the functions used ("*All functions*"), the analysis in the MESA dataset ("*MESA analysis*") and the analysis in cohort 2 ("*Cohort2 Analysis*"). 

There are a few externally developed packages in MATLAB (which do not require installation) in the "*All functions*" folder:

* [ShadedErrorBar](https://uk.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar)
* [ViolinPlot](https://github.com/bastibe/Violinplot-Matlab)
* [DistinguishableColors](https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors#functions_tab)
* For PSG (EDF format) [Reading](https://www.edfplus.info/downloads/index.html)
* Circular statistics toolbox [CircStat](https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)

The other subfolders within contain our self-developed methods described in the manuscript.

> The codes are commented within; If you have any questions about using, or about detailed variable explanations please get in touch with us.

----
#### Analysis codes of the MESA dataset

The analysis scripts are introduced in a sequential order below:

* ***RawPSG_Loading***: This folder contains the MATLAB codes for reading the raw PSG files (*Data_load.m*), as well as participant exclusion checks (*Demo\_Sleep_Checks.m*), and demographic analysis (*Demo_Checks.m*). This step will return a list of Matlab data files for feature extraction in the second step.

* ***Feature Extractions***: This folder contains the code for feature extraction (*fteval_allsbj.m*), using the data coming from the first step.

* ***FeatureSpaceAnalysis***: This folder contains the main analysis code to generate results for Figure 2-3 in the manuscript.
	1. *Preparation*: This folder contains the script (*FeatureDynam_Prep.m*) to prepare the features from the second step towards further group-level feature-space analysis. They include artefact checks on EEG epochs, as well as global z-score normalisation.
	2. *MainVariables\_Group_Figure2*: Group-level analysis of the feature space variables including:
		* Sleep distance `s(t)`: *_SleepDistanceS.m_*
		* State velocity `v(t)`: *_StateVelocityVt.m_*
		* Sleep velocity `v_s(t)`: *_SleepVelocityVs.m_*
		* Bifurcation function fitting of the sleep distance: *_BifFitting\_SleepDistance.m_*
	
		These codes generate the results for parts of Figures 2 and 3.
	
	3. *FeatureAnalyis_Figure3*: Group-level analysis of individual feature dynamics:
	   * FPCA analysis: *_FPCA\_Features.py_*
	   * Granger Causality analysis: *_GrangerCausalityTest.m_*
	   
	   These codes generate the results for part of Figure 3. 
	4. *IndividualAnalysis_Figure2*: Analysis of the individual-level sleep-distance dynamics (for results in Figure 2F-I):
		* *TSPrep_EachIndividual.m*: this code generates the Sleep distance s(t) dynamics for each individual (With some further exclusions);
		* Early warning signal for Critical slowing down analysis: *_EWS\_Evaluations.m_* and *_EWSAnalysis\_Stats.m_*.
		* Bifurcation function fitting on individual-level dynamics analysis: *_Individual\_BifFitting.m_* and *_Individual\_BifFitting\_Analysis.m_*

To see the results, go to the figures folder and visualise them.


-----
#### Analysis for Cohort 2

The analysis scripts are introduced in sequential order below; Analysis of cohort 2 generates results for Figure 4. 

1. xxx

2. xxx

3. xxx

4. *PostProcessing*: The folder post-process all results generated previously for displaying in figures, as well as statistics summary.
	
	* Computing the tipping points from the post-hoc bifurcation function fits for the Sleep distance dynamics from all nights (*_TippingEval\_RawS\_R2Summary.m_*)
	* Adding post-hoc bifurcation function fit tipping point times toward the prediction results (*_AddPostHoc\_ToPrediction.m_*)
	* Summary analysis of the prediction results of the sleep distance dynamics (compared to ground truth); The accuracy was evaluated using the cosine similarity scores; (*_Prediction\_SDynamic.m_*)
	* Summary analysis of the tipping point prediction results (*_PredictionTipping\_Summary.m_*)

To see the results, go to the figures folder and visualise them.

	
















