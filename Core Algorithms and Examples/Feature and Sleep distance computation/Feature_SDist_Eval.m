%% Sleep distance computation and bifurcation model fitting
% This script demonstrate EEG feature computed from a single subject and
% transformed into sleep distance variable (S-variable) and further apply
% bifurcation model fitting

% ** Add all files and subfolder under PSG_example before running the demo
% ** One could skip the feature computation by directly starting from line 91 to load the pre-calculated features and Sleep distance variable 
% ** The feature computation is going to take >10 min in normal laptop case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Important Notice

% This code is simply a demonstration for key aspects of the computational
% pipeline, where some steps (channel exclusions, ECG artefact cleaning)
% are not included:

% Key aspects included:
% 1. Feature computation (with three added features)
% 2. Sleep distance calculation

% Notice 2:
% The feature normalisation (z-score) here utilised the knowledge of
% group-level mean and std of each feature across all participants; For
% Cohort 2 analysis (prediction), this was normalised within individual,
% using values from the training night only; No knowledge from testing
% nights were used as well.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Authors: Junheng Li, Tianyu Wei

%% Load the EEG data and group level feature mean/std

clear all
clc
load EEGProc.mat

% Notice that the EEG has been processed as described in the manuscript,
% and we only used one EEG channel to shorten the runtime as an example
% here. We also cut the EEG to the period 30 min before sleep onset to 10
% minutes after.

%% Parameters
% Define common parameters
Fs = 256;                           % sampling frequency
overlap_min = 0.05;                 % overlap between EEG epochs in min
tper_now = [-30,10];      % periods taken around sleep onset (first two consecutive N2)
overlap_smp = overlap_min*60*Fs;   % number of datapoint for overlap samples (3 s)
epc_size_smp = 2*overlap_smp;       % number of datapoint for a complete epoch (6 s)

%% Divide EEG into epochs and run feature evaluation
% Note that this part is going to take >10 minutes approximately

ft_mat = {};

for i = 1:floor(length(eeg1)/overlap_smp)-1

    epc_idx = overlap_smp*(i-1)+1:overlap_smp*(i-1)+epc_size_smp; % select the epoch to compute

    ft_values_all_ch = []; % Feature from all channels

    % Only compute one channel

    EEG_ch_epc = eeg1(1,epc_idx);
    [ft, ft_dscrp] = sleep_ft_extract_cont(EEG_ch_epc, Fs); % compute the EEG features

    ftraw = cell2mat(struct2cell(ft));

    % Compute additional 3 features
    [ftnew, ft_dscrp_add] = ftadd(EEG_ch_epc, Fs);
    ftnew = cell2mat(struct2cell(ftnew));

    ft_dscrp_all  = [ft_dscrp,ft_dscrp_add];

    ft_values_all_ch = [ftraw;ftnew]';

    ft_mat{end+1} = mean(ft_values_all_ch,1,"omitmissing");

end

%% Sleep distance variable computation 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you want to skip the feature calculation stage (which is long)
% load Features.mat

load GrpFt_MeanStdnow.mat

overlap_min = 0.05;
onset_mediantime = 10;
onsetepcs = onset_mediantime./overlap_min -1;   % Calculate number of median epochs

% concatonate all epoch together
ft_mat_all = vertcat(ft_mat{:});

% Z-score normalize all features to the group level mean and std
ft_mat_all_norm = (ft_mat_all-mean_allft)./std_allft; 

% compute the centroid of the sleep by taking the median of 10 min post onset
ft_onset_centroid = median(ft_mat_all_norm(end-onsetepcs:end,:),1,"omitmissing"); 

% compute the sleep distance by taking the n-dimensional euclidean distance
S_var = pdist2(ft_mat_all_norm,ft_onset_centroid,"euclidean"); 

% smoothing the data 
S_var_smooth = smoothdata(S_var,1,"movmedian",[20,0]); 

xx_smooth = S_var_smooth-min(S_var_smooth);
xx_smooth = xx_smooth';   % Ensure row vector
tvec_plot = -29.9:overlap_min:10;   % Notice causality here, time 0 preceds the first stage N2

figure
plot(tvec_plot,xx_smooth,'LineWidth',2)
box off
xlabel('Time (min)')
ylabel('Sleep distance')
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

%% Save this for the next example, bifurcation fitting

save('SleepDistance_Example.mat','xx_smooth','tvec_plot') % 

