function [base_per,  Onset_per, ifgood, stp_point ,stp_ori ,msg] = AASM_scoring_identify_artnan(scoring, data, ep_size,stage,cont_epc)

% This function is used to identify scoring based on AASM criterion of EEG
% data, and extractes first minute N1 recording, and first minute N2 recording
%
% Author: Junheng Li
%
% Input:
% scoring: the scoring of epochs
% data: relevant EEG recording, can be a NxM matrix where N indicate
% channels and M indicate number of data points
% ep_size: epoch size used in scoring (number of data points but not
% second)
% stage: string indicating whether to extract stage N1 onset or N2 onset
% cont_epc: Number of continous same stage epochs to extract 
%
% Output:
% base_per: awake state recording extracted with specified length
% Onset: specified stage state recording extracted with specified length
% ifgood: A logical value telling if the data is good to use, if there is
% not enough length of recording it is 0, and the data will be discarded
% stp_point: the start point of each period if the data is good
% msg: information of scoring (why not good)
%% Log of code
% All previous changes: modified the code to consider artefact epochs
% (marked as NaNs). Also added outputs of two different types of sleep
% onset markers (original location with NaNs and NaN-removed locations)

% 30 July 2021, changing the maximum loop size to fit the codes;

% 03 August 2021, changing again the maximum loop size in order to
% eliminate maximum index errors

%%

base_per = [];
Onset_per = [];
ifgood = 1;
stp_point = [];
stp_ori = [];

max_num_epc = floor(size(data,2)/ep_size);
scoring_clean = scoring;
scoring_clean(isnan(scoring)) = [];
% scoring_clean = scoring_clean(1:max_num_epc);    % Limit the size of scoring to that of the signal

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find continuous awake baseline starts
base_start = find(scoring_clean == 0, 1);
if isempty(base_start)    % No such stage
    ifgood = 0;
    msg = 'No awake stage found';
    return
end
id_Onset = 1;
max_loop = length(find(scoring_clean == 0))-1;
if base_start >= max_num_epc-cont_epc+2
    msg = 'No awake stage found';
    ifgood = 0;
    return
end
while sum(scoring_clean(base_start:base_start+cont_epc-1)) ~= 0 && id_Onset<max_loop
    id_Onset = id_Onset+1;
    base_start = find(scoring_clean == 0, id_Onset);
    if length(base_start)<id_Onset
        ifgood = 0;
        msg = 'No continuous 1min awake period found';
        return
    end
    base_start = base_start(id_Onset);
    if base_start >= max_num_epc-cont_epc+2
        ifgood = 0;
        msg = 'No continuous 1min awake period found';
        return
    end
end

idx_noart = find((isnan(scoring)) == 0);
base_real = idx_noart(base_start);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the continous stage start

% Make sure that the N2 period taken are behind the awake period taken
if base_start+cont_epc >= max_num_epc-cont_epc+2
    ifgood = 0;
    msg = 'Awake found at the end of recording';
    return
end
scoring_clean = scoring_clean(base_start+cont_epc:end);

if stage == 'N1'
    stage_num = 1;
    sum_cont = stage_num*cont_epc;
elseif stage == 'N2'
    stage_num = 2;
    sum_cont = stage_num*cont_epc;
end

Stage_start = find(scoring_clean == stage_num, 1);
if isempty(Stage_start)    % No such stage
    ifgood = 0;
    msg = 'No sleep stage after continuous awake';
    return
end
id_Onset = 1;
max_loop = length(find(scoring_clean == stage_num))-1;
if Stage_start >= max_num_epc-cont_epc+2
    ifgood = 0;
    msg = 'No sleep stage after continuous awake';
    return
end
% Find first continuous stage minute
while (sum(scoring_clean(Stage_start:Stage_start+cont_epc-1)) < sum_cont || sum(scoring_clean(Stage_start:Stage_start+cont_epc-1))>6)&& id_Onset<max_loop
    id_Onset = id_Onset+1;
    Stage_start = find(scoring_clean == stage_num, id_Onset);
    if length(Stage_start)<id_Onset
        ifgood = 0;
        msg = 'No continuous 1min sleep stage found';
        return
    end
    Stage_start = Stage_start(id_Onset);
    if Stage_start >= max_num_epc-cont_epc+2
        ifgood = 0;
        msg = 'No continuous 1min sleep stage found';
        return
    end
end
if sum(scoring_clean(Stage_start:Stage_start+cont_epc-1)) < sum_cont || sum(scoring_clean(Stage_start:Stage_start+cont_epc-1))>6
    ifgood = 0;
    msg = 'No continuous 1min sleep stage found';
    return
end
% Correct the right index
Stage_start = Stage_start+base_start+cont_epc-1;
stage_real = idx_noart(Stage_start);

%% Take the period

base_per = data(:,(base_start-1)*ep_size+1:(base_start+cont_epc-1)*ep_size);
Onset_per = data(:,(Stage_start-1)*ep_size+1:(Stage_start+cont_epc-1)*ep_size);
stp_point = [base_real,stage_real];
stp_ori = [base_start,Stage_start];
msg = 'All good';

