%% Demographic info corrections for sleep-onset latency to keep consistent

%% Preparation

clear all;
clc

%% General parameters

cont_def = 2*30;    % Length of time to define a continuous stage

quality_th = 4;    % Minimum acceptable EEG quality

%% Loading data

for part = 1:52

load(['data_all_mesa_part',num2str(part),'.mat'])

%% Parameters

num_sbjs = length(data_all_mesa);    % Number of subjects


%% Data processing & epoch extraction

% Using sham data as a first step
Stage = 'N2';    % Stage onset to study

base_epcs_all = [];
Stage_epcs_all = [];

labels_base = {};
labels_Stage = {};

discard_info = zeros(num_sbjs,3);
discard_reason = cell(num_sbjs,1);
original_ids = zeros(num_sbjs,1);

sleep_info = zeros(num_sbjs,2);
sleep_info_artremoved = zeros(num_sbjs,2);
sleep_info_oriscore_nonan = NaN(num_sbjs,2);
time_to_sleep = zeros(num_sbjs,1);
time_to_sleep_correct = NaN(num_sbjs,1);
ecgtrace = cell(num_sbjs,1);
diff_list = [];

for sbj = 1:num_sbjs
   
    score_size = data_all_mesa{sbj}.scoreEpoch;   % Scoring epoch size in s
    Fs = data_all_mesa{sbj}.EEG_sampling_freq;   % Sampling rate
    
    channels = data_all_mesa{sbj}.EEG_channel;    % Channels of EEG
    
    original_ids(sbj) = data_all_mesa{sbj}.Orignal_id;
   
        scoring_now = data_all_mesa{sbj}.scoringLabels;    % Current scoring
        
        eeg_mat_now = data_all_mesa{sbj}.EEG_mat;
        ecg_now =  data_all_mesa{sbj}.ECG_sig{1}';
        
        % Remove the start unstable signal part
%         N_rm = floor(start_remove/score_size);
%         [eeg_mat_rm,score_rm] = remove_start(N_rm,eeg_mat_now,scoring_now,score_size*Fs);
        
        N_cont = ceil(cont_def/score_size);
        
        % Check EEG quality
        EEG_quality_now =  data_all_mesa{sbj}.EEG_quality;
        num_chs = length(channels);
        delete_ch = zeros(num_chs,1);
        for ch = 1:num_chs
            if EEG_quality_now(ch) < quality_th
                disp(['Subject ',num2str(data_all_mesa{sbj}.Orignal_id), ' channel ',channels{ch}, ' discarded due to poor quality'])
                discard_info(sbj,ch) = 1;
                discard_reason{sbj,1} = 'poor quality';
                delete_ch(ch) = 1;
                continue;
            end
        end
        idx_rm_ch = (delete_ch == 0);
        eeg_mat_now = eeg_mat_now(idx_rm_ch,:);
        
        % Check if all channels are removed, then skip this subject
        if isempty(eeg_mat_now)
            continue
        end

        % General artefact rejection based on scoring epochs
        [eeg_mat_artnan,ecg_artnan,~,~, score_artfree, score_clean] = genart_rm_artnan(eeg_mat_now,ecg_now,Fs,score_size,2,scoring_now);
        score_ori = scoring_now(1:length(score_artfree));

        % ECG artefact removal
        [eeg_mat_clean,ori_trace,res_ecg] = ecgart_rm(eeg_mat_artnan,ecg_artnan,Fs,0.02);
        ecgtrace{sbj}.original = ori_trace;
        ecgtrace{sbj}.residual = res_ecg;
        
        % Scoring identification
        [base_per_mat,~,ifgood1,stp_points,nnn,msg] = AASM_scoring_identify_artnan(score_artfree,eeg_mat_clean,score_size*Fs,Stage,N_cont);
        [~, Stage_per_mat,ifgood2,stp_points_ori,nnn_ori,msg] = AASM_scoring_identify_artnan(score_ori,eeg_mat_clean,score_size*Fs,Stage,N_cont);
        if ~ifgood1      % No such stage onset
            disp(['Subject ',num2str(data_all_mesa{sbj}.Orignal_id), ' ', ' is discarded due to no useful scoring epochs'])
            discard_info(sbj,:) = 1;
            discard_reason{sbj,1} = 'Useless scoring';
            continue;
        end
        if ~ifgood2      % No such stage onset
            disp(['Subject ',num2str(data_all_mesa{sbj}.Orignal_id), ' ', ' is discarded due to no useful scoring epochs'])
            discard_info(sbj,:) = 1;
            discard_reason{sbj,1} = 'Useless scoring';
            continue;
        end
        if sum(score_artfree(1:stp_points(1)),'omitnan') > 0     % Remove subjects that did not start from awake
            discard_info(sbj,:) = 1;
            discard_reason{sbj,1} = 'No awake start';
            time_to_sleep(sbj) = 0;
            time_to_sleep_correct(sbj) = 0;
            disp(['Subject ',num2str(data_all_mesa{sbj}.Orignal_id), ' ', ' is discarded due to not starting from awake'])
            continue
        end
        sleep_info(sbj,:) = stp_points;
        sleep_info_artremoved(sbj,:) = nnn;
        sleep_info_oriscore_nonan(sbj,:) = stp_points_ori;
        time_to_sleep(sbj) = diff(stp_points)*score_size;   % Time to sleep in seconds
        time_to_sleep_correct(sbj) = (stp_points_ori(2) - stp_points(1))*score_size;
        
        % New checks
        if stp_points_ori(2)<= stp_points(1)           % Contradictory definition
            disp(['Subject ',num2str(data_all_mesa{sbj}.Orignal_id), ' ', ' inconsistent bedtime and sleep-onset after redefinition'])
            time_to_sleep_correct(sbj) = NaN;
        end
        if stp_points_ori(2) ~= stp_points(2)          % Different SO point check due to NaNs
            diff_list = [diff_list,sbj];
        end


        
        if time_to_sleep(sbj)>5400 || time_to_sleep(sbj)<180   % Remove too long fall asleep time
            discard_info(sbj,:) = 1;
            discard_reason{sbj,1} = 'Strange asleep time';
            disp(['Subject ',num2str(data_all_mesa{sbj}.Orignal_id), ' ', ' is discarded due to too long asleep time'])
            continue
        end
       
       
    
        disp(['Subject ',num2str(data_all_mesa{sbj}.Orignal_id),' Finished' ])
end
cd('/home/junheng/Dropbox (UK Dementia Research Institute)/Data proc and analysis/Analysis MESA/General demographic/Subject Demo Correction July2023')
save(['Subject_info','_part_',num2str(part),'.mat'],'discard_info','original_ids','ecgtrace','sleep_info','sleep_info_artremoved','time_to_sleep','discard_reason','diff_list','sleep_info_oriscore_nonan','time_to_sleep_correct')



end