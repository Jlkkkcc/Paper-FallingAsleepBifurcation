%% MESA Raw PSG data loadings

% Author: Junheng Li

% This code is written for extracting NSRR data from Original dataset and
% trying to reorganinze them;
% This code also extracts the .XML file for signal annotation, and aligns
% it with the data

% Demographic information should be converted to .mat manually before this
% step for "MESA_DemoInfo_raw.mat" / Manual checking for PSG files to ensure the raw subject indexes
% were loaded for "MESA_EDFid_raw.mat";

%% Preparation

clear all
clc

load MESA_EDFid_raw.mat
load MESA_DemoInfo_raw.mat
load MESA_ECGquality_raw.mat
load MESA_EEGquality_raw.mat

num_sbjs = length(MESA_EDFid);     % Number of subjects available in total

%% File detection (remove non-existing subjects

fdr = '/mnt/data1TB/mesa/polysomnography/edfs';
cd(fdr)

idx_rawselect = zeros(num_sbjs,1);

% Detect non-existing EDFs
for sbj = 1:num_sbjs
    prefix = 'mesa-sleep-';
    postfix = '.edf';
    
    % File name generator
    if MESA_EDFid(sbj)<10
        fid = [prefix,'000',num2str(MESA_EDFid(sbj)),postfix];
    elseif MESA_EDFid(sbj)<100 && MESA_EDFid(sbj)>=10
        fid = [prefix,'00',num2str(MESA_EDFid(sbj)),postfix];
    elseif MESA_EDFid(sbj)<1000 && MESA_EDFid(sbj)>=100
        fid = [prefix,'0',num2str(MESA_EDFid(sbj)),postfix];
    else
        fid = [prefix,num2str(MESA_EDFid(sbj)),postfix];
    end
    
    fgood = fopen(fid);
    if fgood == -1
        MESA_EDFid(sbj) = NaN;
        Gender_raw(sbj) = NaN;
        Ages_raw(sbj) = NaN;
        Race_raw(sbj) = NaN;
        ECG_quality_raw(sbj) = NaN;
        EEG_quality_C4_raw(sbj) = NaN;
        EEG_quality_Cz_raw(sbj) = NaN;
        EEG_quality_Fz_raw(sbj) = NaN;
        continue
    end
    idx_rawselect(sbj) = 1;
end

MESA_EDFid(isnan(MESA_EDFid)) = [];   % Remove NaNs

Gender_raw(isnan(Gender_raw)) = [];
Gender = Gender_raw;

Ages_raw(isnan(Ages_raw)) = [];
Age = Ages_raw;

Race_raw(isnan(Race_raw)) = [];
Race = Race_raw;

ECG_quality_raw(isnan(ECG_quality_raw)) = [];
ECG_quality = ECG_quality_raw;

EEG_quality_C4_raw(isnan(EEG_quality_C4_raw)) = [];
EEG_quality_C4 = EEG_quality_C4_raw;

EEG_quality_Cz_raw(isnan(EEG_quality_Cz_raw)) = [];
EEG_quality_Cz = EEG_quality_Cz_raw;

EEG_quality_Fz_raw(isnan(EEG_quality_Fz_raw)) = [];
EEG_quality_Fz = EEG_quality_Fz_raw;

num_sbjs = length(MESA_EDFid);

save("MESA_EDFid_exist.mat","MESA_EDFid");
save("MESA_Demoinfo.mat","Age","Gender","Race");
save("MESA_ECGquality.mat","ECG_quality");
save("MESA_EEGquality.mat","EEG_quality_Fz","EEG_quality_C4","EEG_quality_Cz")


%% EDF data loading
% The entire dataset is too large, therefore divide it into parts

clear all
clc

load MESA_EDFid_exist.mat
num_sbjs = length(MESA_EDFid);     % Number of subjects available in total

% Load demographic information
load MESA_DemoInfo.mat
load MESA_ECGquality.mat
load MESA_EEGquality.mat

% Part x
max_sbj = 40;

for part = 1

EEG_header_strs = {'EEG1','EEG2','EEG3'};    % The lables for the 3 EEG channels
EEG_chs = {'Fz-Cz','Cz-Oz','C4-M1'};

% Fact of sampling frequency
samp_freq_EEG = 256;
samp_freq_ECG = 256;
hours_used = 6;
time_to_take = hours_used*60*60;    % Time period of signal to take
datapoints_used = time_to_take*samp_freq_EEG;


ECG_header_strs = {'EKG'};               % ECG channel label

fdr = '/mnt/data1TB/mesa/polysomnography/edfs';
cd(fdr)

st_idx = (part-1)*max_sbj;
end_idx = part*max_sbj;
if end_idx > num_sbjs
    max_sbj = num_sbjs - st_idx;
end

data_all_mesa = cell(max_sbj,1);

for sbj = 1:max_sbj
    prefix = 'mesa-sleep-';
    postfix = '.edf';
    
    % File name generator
    if MESA_EDFid(st_idx+sbj)<10
        fid = [prefix,'000',num2str(MESA_EDFid(st_idx+sbj)),postfix];
    elseif MESA_EDFid(st_idx+sbj)<100 && MESA_EDFid(st_idx+sbj)>=10
        fid = [prefix,'00',num2str(MESA_EDFid(st_idx+sbj)),postfix];
    elseif MESA_EDFid(st_idx+sbj)<1000 && MESA_EDFid(st_idx+sbj)>=100
        fid = [prefix,'0',num2str(MESA_EDFid(st_idx+sbj)),postfix];
    else
        fid = [prefix,num2str(MESA_EDFid(st_idx+sbj)),postfix];
    end
    data_all_mesa{sbj}.Orignal_id = MESA_EDFid(st_idx+sbj);
    
    % Reading corresponding edf of current subject
    [data_all_mesa{sbj}.Header, signalHeader, signalCell] = blockEdfLoad(fid);
      
    % ECG cell location
    len_header = length(signalHeader);
    kk = 1;
    num_ECG = length(ECG_header_strs);
    ECG_idx = zeros(num_ECG,1);
    while kk < num_ECG+1
        for i = 1:len_header
            sig_lab = signalHeader(i).signal_labels;
            if strcmp(ECG_header_strs{kk},sig_lab)
                ECG_idx(kk) = i;
                break
            end
        end
        kk = kk+1;
    end
    
    % EEG cell location
    kk = 1;
    num_EEG = length(EEG_header_strs);
    EEG_idx = zeros(num_EEG,1);
    while kk < num_EEG+1
        for i = 1:len_header
            sig_lab = signalHeader(i).signal_labels;
            if strcmp(EEG_header_strs{kk},sig_lab)
                EEG_idx(kk) = i;
                break
            end
        end
        kk = kk+1;
    end
    
    % EEG signal storing
    len_EEGsig = datapoints_used;
    length(signalCell{EEG_idx(1)})
    if length(signalCell{EEG_idx(1)})<datapoints_used
        len_EEGsig = length(signalCell{EEG_idx(1)});
    end
    EEG_mat = zeros(num_EEG,len_EEGsig);
    for i = 1:num_EEG
        data_all_mesa{sbj}.EEG_sig{i} = signalCell{EEG_idx(i)}(1:len_EEGsig);
        EEG_mat(i,:) = signalCell{EEG_idx(i)}(1:len_EEGsig);
        data_all_mesa{sbj}.EEG_info{i} = signalHeader(EEG_idx(i));
    end
    data_all_mesa{sbj}.EEG_mat = EEG_mat;   % Matrix form EEG data
    data_all_mesa{sbj}.EEG_channel = EEG_chs;    % Corresponding EEG channels
    data_all_mesa{sbj}.EEG_quality(1) = EEG_quality_Fz(st_idx+sbj);
    data_all_mesa{sbj}.EEG_quality(2) = EEG_quality_Cz(st_idx+sbj);
    data_all_mesa{sbj}.EEG_quality(3) = EEG_quality_C4(st_idx+sbj);
    data_all_mesa{sbj}.EEG_sampling_freq = samp_freq_EEG;
    
    % ECG signal storing
    for i = 1:num_ECG
        data_all_mesa{sbj}.ECG_sig{i} = signalCell{ECG_idx(i)}(1:len_EEGsig);
        data_all_mesa{sbj}.ECG_info{i} = signalHeader(i);
        data_all_mesa{sbj}.ECG_quality(i) = ECG_quality(st_idx+sbj);
    end
    data_all_mesa{sbj}.ECG_sampling_freq = samp_freq_ECG;
    
    % Demographic information
    
    data_all_mesa{sbj}.age = Age(st_idx+sbj);
    data_all_mesa{sbj}.race = Race(st_idx+sbj);
    data_all_mesa{sbj}.gender = Gender(st_idx+sbj);
    
    disp(['subject No.',num2str(MESA_EDFid(st_idx+sbj)),' is loaded'])

end

fdr = '/mnt/data1TB/mesa/polysomnography/annotations-events-nsrr';
cd(fdr)
% Reading and aligning annotations
for sbj = 1:max_sbj
    prefix = 'mesa-sleep-';
    postfix = '-nsrr.xml';
    
     % File name generator
    if MESA_EDFid(st_idx+sbj)<10
        fid = [prefix,'000',num2str(MESA_EDFid(st_idx+sbj)),postfix];
    elseif MESA_EDFid(st_idx+sbj)<100 && MESA_EDFid(st_idx+sbj)>=10
        fid = [prefix,'00',num2str(MESA_EDFid(st_idx+sbj)),postfix];
    elseif MESA_EDFid(st_idx+sbj)<1000 && MESA_EDFid(st_idx+sbj)>=100
        fid = [prefix,'0',num2str(MESA_EDFid(st_idx+sbj)),postfix];
    else
        fid = [prefix,num2str(MESA_EDFid(st_idx+sbj)),postfix];
    end
    
    ifannot = fopen(fid);
    if ifannot == -1
        data_all_mesa{sbj}.annot = -1;
        disp(['No relevant annotation found for subject No.',num2str(MESA_EDFid(st_idx+sbj))])
        continue
    else
        data_all_mesa{sbj}.annot = 1;
        obj = loadPSGAnnotationClass(fid);
        obj_annot = obj.loadFile;           %Extract the contents
        
        % Extract Arousal information
        data_all_mesa{sbj}.ArousalsStart = cell2mat(obj_annot.ArousalsStart);
        data_all_mesa{sbj}.ArousalsDuration = cell2mat(obj_annot.ArousalsDuration);
        
        % Pre-process Sleep stages
        stages_now = cell2mat(obj_annot.SleepStages);
        score_epoch_len = obj_annot.EpochLength;    % Scoring length
        stages_dur = cell2mat(obj_annot.SleepStagesDuration);
        stages_start = cell2mat(obj_annot.SleepStagesStart);
        
        len_score = length(stages_start);
        new_score_vec = [];    % New scoring vector for each epoch
        new_score_start = [];   % Starting point of each scoring unit
        num_epc_curscore = zeros(len_score,1);
        for i = 1:len_score
            if i <len_score
                if stages_start(i)+stages_dur(i) == stages_start(i+1)         % Continuous
                    num_epc_curscore(i) = floor(stages_dur(i)/score_epoch_len);
                    new_score_vec = [new_score_vec,ones(1,num_epc_curscore(i))*str2double(stages_now(i))];
                else
                    num_empty_epc = (stages_start(i+1) - (stages_start(i)+stages_dur(i)))/score_epoch_len;
                    num_epc_curscore(i) = floor(stages_dur(i)/score_epoch_len);
                    new_score_vec = [new_score_vec,ones(1,num_epc_curscore(i))*str2double(stages_now(i))];
                    new_score_vec = [new_score_vec,nan(1,num_empty_epc)];
                end
            else
                num_epc_curscore(i) = floor(stages_dur(i)/score_epoch_len);
                new_score_vec = [new_score_vec,ones(1,num_epc_curscore(i))*str2double(stages_now(i))];
            end
        end
        score_start = stages_start(1);
        score_len_total = length(new_score_vec)*score_epoch_len;
        
        % Storing necessary information
        data_all_mesa{sbj}.scoringLabels = new_score_vec;
        data_all_mesa{sbj}.scoringStart = score_start;
        data_all_mesa{sbj}.scoringPeriod = score_len_total;
        data_all_mesa{sbj}.scoreEpoch = score_epoch_len;
        
    end
    
    disp(['subject No.',num2str(MESA_EDFid(st_idx+sbj)),' is annotated'])
    
end

fdr = '/mnt/data1TB/Data proc and analysis/Data loading MESA/Raw PSG data all';
cd(fdr)
save(['data_all_mesa_part',num2str(part),'.mat'],'data_all_mesa')

disp(['Part No.',num2str(part),' is successfully saved'])
clear data_all_mesa 

end











