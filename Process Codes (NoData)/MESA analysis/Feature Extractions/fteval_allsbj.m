%% Feature extraction continuous
%
% Author: Junheng LI
% 
% This code is written for continuous feature extraction on the MESA
% dataset;

%% Preparations

% Subject group
clear all
clc
load Sbjs_filtered_time.mat

sbjs_ids = sbj_remained;

%% Continuous feature evaluation pipelines

% Original analysis parameters
Stage = 'N2';    % Stage onset to study
num_per_part = 40;
part_num_old = 1;

cont_def = 60;
quality_th = 4;
num_sbjs = length(sbjs_ids);

% Continuous extraction parameters

Fs = 256;   % Sampling rate

epc_len = 6;            % Epoch length for feature evaluation
ovlap = 0.5;              % Overlap between continous epochs
post_sleep_per = 20;    % Post sleep period for feature evaluation

discard_idx = [];
num_ftepcs_all = NaN(num_sbjs,1);
stp_ftepcs_all = [];
epcs_to_asleep = NaN(num_sbjs,1);
oriscore_asleep = [];

% Evaluation and pipelines

for i = 1:num_sbjs
    
    part_num = floor(sbjs_ids(i)/(num_per_part+0.00001))+1;
    if i ~=1
        if part_num ~= part_num_old
            load(['data_all_mesa_part',num2str(part_num),'.mat'])
            part_num_old = part_num;
        end
    elseif i == 1
        load(['data_all_mesa_part',num2str(part_num),'.mat'])
    end
    
    sbj = sbjs_ids(i)- (part_num-1) * num_per_part;
    
    score_size = data_all_mesa{sbj}.scoreEpoch;   % Scoring epoch size in s
    Fs = data_all_mesa{sbj}.EEG_sampling_freq;   % Sampling rate
    
    channels = data_all_mesa{sbj}.EEG_channel;    % Channels of EEG

    scoring_now = data_all_mesa{sbj}.scoringLabels;    % Current scoring
    
    eeg_mat_now = data_all_mesa{sbj}.EEG_mat;
    N_cont = ceil(cont_def/score_size);
    ecg_now =  data_all_mesa{sbj}.ECG_sig{1}';
    
    % Check EEG quality
    EEG_quality_now =  data_all_mesa{sbj}.EEG_quality;
    num_chs = length(channels);
    delete_ch = zeros(num_chs,1);
    for ch = 1:num_chs
        if EEG_quality_now(ch) < quality_th
            disp(['Subject ',num2str(data_all_mesa{sbj}.Orignal_id), ' channel ',channels{ch}, ' discarded due to poor quality'])
            
            delete_ch(ch) = 1;
            continue;
        end
    end
    idx_rm_ch = (delete_ch == 0);
    eeg_mat_now = eeg_mat_now(idx_rm_ch,:);
    
    % Check if all channels are removed, then skip this subject
    if isempty(eeg_mat_now)
        discard_idx = [discard_idx;i];
        continue
    end
    len_score = floor(length(ecg_now) / (score_size*Fs));
    score_ori = scoring_now(1:len_score);     % Original scoring without NaN artefact marks

    [eeg_mat_now,ori_trace,res_ecg] = ecgart_rm(eeg_mat_now,ecg_now,Fs,0.02);
    
    % Scoring check
    [~, ~,ifgood,stp_points,stp_points_ori,~] = AASM_scoring_identify_artnan(score_ori,eeg_mat_now,score_size*Fs,Stage,N_cont);
    if ~ifgood      % No such stage onset
        disp(['Subject ',num2str(data_all_mesa{sbj}.Orignal_id), ' ', ' is discarded due to no useful scoring epochs'])
        continue;
    end
    
    % Specific selection criterion
    slepc_toasleep = diff(stp_points_ori)+1;      % Number of scoring epochs to asleep
    % Detect start of identification
    if stp_points_ori(1) ~=1
        disp(['Sleep scoring not starting from beginning, Subject No.',num2str(sbjs_ids(i))])
        discard_idx = [discard_idx;i];
        continue
    end
    epcs_to_asleep(i) = slepc_toasleep;
    score_to_asleep = score_ori(stp_points_ori(1):stp_points_ori(2)+post_sleep_per);
    oriscore_asleep{i} = score_to_asleep;

    % Feature evaluation
    idx_realch = 0;
    for ch = 1:length(channels)
        
        if delete_ch(ch)   % Check if this channel is deleted
            continue
        end
        idx_realch = idx_realch+1;
        
        eeg_now = eeg_mat_now(idx_realch,:);
        
        start_p = (stp_points_ori(1)-1)*score_size*Fs+1;
        stop_p = (stp_points_ori(2)+post_sleep_per)*score_size*Fs ;
        eeg_ck = eeg_now(start_p:stop_p);      % Entire period from start to the post-asleep period
        eeg_asleepper_all{i}.ch{ch} = eeg_ck;
        
        len_now = length(eeg_ck);
        epc_len_data = epc_len*Fs;
        ovlap_data = floor(epc_len_data*ovlap);
        rm_data = epc_len_data-ovlap_data;
        
        num_epcs_total = floor((len_now - ovlap_data)/rm_data);        % Total number of epochs given overlapping
        num_ftepcs_all(i) = num_epcs_total;
        
        stp_ftepcs = zeros(num_epcs_total,1);
        
        for jj = 1:num_epcs_total
            stp = (jj-1)*rm_data +1;            % Starting timepoint of current epoch
            stp_ftepcs(jj) = stp;
            data_unit = eeg_ck(stp:stp+epc_len_data-1);     % Current data chunk
            if max(data_unit) == 0.5    % Data exceding the Amplifier accuracy
                ft_sbj{i}.ck{jj}.ft_ch(ch,:) = NaN(1,47);
            else
                [ftt,ft_dscrp] = sleep_ft_extract_cont(data_unit,Fs);
                
                ftt1 = struct2cell(ftt);
                
                ft_sbj{i}.ck{jj}.ft_ch(ch,:) = transpose(cell2mat(ftt1));
            end
            
        end

    end
    stp_ftepcs_all{i} = stp_ftepcs;
    
    disp(['Subject No.',num2str(sbjs_ids(i)),' successfully calculated'])


end

clear data_all_mesa
save('fteval_EWSallsbj_20Dec.mat')






