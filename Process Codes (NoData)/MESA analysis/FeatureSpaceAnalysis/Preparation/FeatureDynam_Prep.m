%% Preparation

% This code runs before feature space analysis to set up for bootstrapping
% & artefact controls for group-level analysis

% This code also returns the group level dynamics for each feature (47
% in total)


%% Data loading
clear all
clc
load 'fteval_EWSallsbj_20Dec.mat'

%% Artefact check based on shorter epochs and EEG extracted;
%
% This part has been modified from previous continuous feature dynamical
% analysis, to increase the feature sampling rate, and for a much better
% dynamical analysis results;
%
% Notice that the EEG data used for artefact check here is the falling
% asleep period + 10 minutes post-asleep period
%
% Revised analysis to use overlapping running window data

% Epoch definition
epc_len = 6;           % Epoch length in seconds
epc_lend = epc_len*Fs;                % Epoch length in data points
ovlap = 0.5;
jmp_len = epc_len - epc_len*ovlap;
jmp_lend = epc_lend - floor(epc_lend*ovlap);

% Artefact specific parameters
nrep = 2;
mul_rms = 2.5;  

sbjs_empty = [];         % A note of empty subjects removed
totalepc_sbj = zeros(num_sbjs,1);

for sbj = 1:num_sbjs

    eeg_sbj_now = eeg_asleepper_all{sbj};
    if isempty(eeg_sbj_now)
        sbjs_empty = [sbjs_empty;sbj];
        continue
    end

    max_ch = length(eeg_sbj_now.ch);          % Current maximum channel number

    % Loop over channel
    isart_all = [];

    for ch = 1:max_ch

        if isempty(eeg_sbj_now.ch{ch})
            continue
        end
        eeg_chnow = eeg_sbj_now.ch{ch};
        len_eeg = length(eeg_chnow);

        num_epcs_max = floor((len_eeg-epc_lend)/(jmp_lend))+1;       % Maximum number of epochs
        isart_ch = zeros(1,num_epcs_max);
        rms_ch = zeros(1,num_epcs_max);
        stps_sbj = zeros(1,num_epcs_max);

        for nep = 1:num_epcs_max
            stps_now = (nep-1)*jmp_lend+1;
            eeg_epc = eeg_chnow(stps_now:stps_now+epc_lend-1);    % The current EEG epoch
            rms_ch(nep) = rms(eeg_epc);
            stps_sbj(nep) = stps_now;
        end
        % Recursive artefact detection threshold
        for j = 1:nrep
            med_rmsch = median(rms_ch,'omitnan');
            isart_ch_now = (rms_ch>(mul_rms*med_rmsch));
            rms_ch(isart_ch_now) = NaN;

        end
        isart_ch = isnan(rms_ch);
        isart_all = [isart_all;isart_ch];

    end
    stps_all{sbj} = stps_sbj;
    isart_sbj = sum(isart_all,1)>0;
    isart_allsbj{sbj} = isart_sbj;
    totalepc_sbj(sbj) = length(isart_sbj);
    
end

%% Feature dynamical time series construction and bootstrap resampling

% The feature values used here are evaluated during the EWS analysis
% The boostrap method is consistent with previous;

% To find the feature values, the indexes used in EWS analysis are matched
% with indexes here.

N_min = 200;         % Minimum number of samples
post_sleep_per = 21;
max_epctotal = max(totalepc_sbj);
epc_postasleep = floor((post_sleep_per*score_size-epc_len)/jmp_len)+1;        % Number of epochs for post-asleep period
ftmark = NaN(num_sbjs,max_epctotal);
ftidx_orid = NaN(num_sbjs,max_epctotal);

for sbj = 1:num_sbjs

    if ismember(sbjs_empty,sbj)
        continue
    end
    
    stps_sbj_new = stps_all{sbj};
    stps_ori = stp_ftepcs_all{sbj};
    % Matching indexes, the EWS version has higher sampling rate
    ftepc_idx_now = find(ismember(stps_ori,stps_sbj_new)==1);
    
    len_ft_now = length(ftepc_idx_now);
    epc_pre = len_ft_now - epc_postasleep;        % Number of epochs prior sleep onset;

    idx_st = max_epctotal - len_ft_now +1;
    ftmark(sbj,idx_st:end) = 1-isart_allsbj{sbj};
    ftidx_orid(sbj,idx_st:end) = ftepc_idx_now;
end

num_samps_time = sum(ftmark,1,'omitnan');
time_start_sampenough = find(num_samps_time>=N_min,1);
disp('Time to sleep onset in minutes:')
((max_epctotal-time_start_sampenough-epc_postasleep-1)*jmp_len+epc_len)/60

% Indexes for bootstraping
samp_enou = (sum(ftmark,1,'omitnan')>=N_min);

for jj = 1:max_epctotal
    
    if jj<time_start_sampenough
        idx_bts{jj} = [];
        continue
    end

    ftidx_now = ftmark(:,jj);
    idxft_all_now = find(ftidx_now==1);
    if length(idxft_all_now)<N_min      % Avoid rare cases after the point the data sample is not enough
        idx_bts{jj} = idxft_all_now;  

    else
        rng(jj)        % Random number seed
        rdx = randperm(length(idxft_all_now),N_min);
        idx_bts{jj} = idxft_all_now(rdx); 
    end

end

%% Feature dynamics plot

num_ft = 50;    % Before revision, 47 features
max_epctotal = max(totalepc_sbj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post 10min norm
onset_mark = (max_epctotal - epc_postasleep+1):(max_epctotal-10);
ftall_mat_allnorm = NaN(num_ft,max_epctotal);

% Feature vector construction
for nft = 1:num_ft
    
    ftmat_this = NaN(num_sbjs,max_epctotal);     % Empty feature matrix for this feature

    for nsbj = 1:num_sbjs
        
        oridxft = ftidx_orid(nsbj,:);     % The feature indexes of previous EWS calculation

        for jj = 1:max_epctotal
            
            if isnan(ftmark(nsbj,jj))            % To remove the artefact and the start empty part
                continue
            end
            ftnow = ft_sbj{nsbj}.ck{oridxft(jj)}.ft_ch(:,nft);    % The current feature values
            ftnow (ftnow == 0) = NaN;
            ftmat_this(nsbj,jj) = mean(ftnow,'omitnan');       % Average features across channels

        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % zscore normalisation globally first
    ftmat_this_zsnorm = (ftmat_this - mean(ftmat_this(:),'omitnan'))./std(ftmat_this(:),'omitnan');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalisation to sleep onset median
    ftonset_this_norm = ftmat_this_zsnorm(:,onset_mark);
    medftonset_this = median(ftonset_this_norm,2,'omitnan');

    ftmat_this_zsnorm_onsetnorm = ftmat_this_zsnorm;
    for nsbj = 1:num_sbjs
        ftmat_this_zsnorm_onsetnorm(nsbj,:) = (ftmat_this_zsnorm(nsbj,:) - medftonset_this(nsbj));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bootstraping
    ftmat_this_zsnorm_onsetnorm_bts = NaN(num_sbjs,max_epctotal);
    for nt = 1:length(idx_bts)
        
        if isempty(idx_bts{nt})
            ftall_mat_allnorm(nft,nt) = NaN;
        else
            ftall_mat_allnorm(nft,nt) = mean(ftmat_this_zsnorm_onsetnorm(idx_bts{nt},nt),'omitnan');
            ftmat_this_zsnorm_onsetnorm_bts(idx_bts{nt},nt) = ftmat_this_zsnorm_onsetnorm(idx_bts{nt},nt);
        end

    end

%     ftall_mat_allnorm(nft,:) = mean(ftmat_this_zsnorm_onsetnorm,1,'omitnan');
    ftall_mat_allsbj{nft} = ftmat_this_zsnorm_onsetnorm;
    ftall_mat_allsbj_bts{nft} = ftmat_this_zsnorm_onsetnorm_bts;

end

%% Compute global mean and std for each feature for normalisation purposes later
% This is time costly so pre-computed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute global mean values and standard deviation for normalisation
% This is computed each feature wise for all subjects;
mean_allft = NaN(1,num_ft);
std_allft = NaN(1,num_ft);
for nft = 1:num_ft
    
    ftmat_this = NaN(num_sbjs,max_epctotal);     % Empty feature matrix for this feature

    for nsbj = 1:num_sbjs
        
        oridxft = ftidx_orid(nsbj,:);     % The feature indexes of previous EWS calculation

        for jj = 1:max_epctotal
            
            if isnan(ftmark(nsbj,jj))            % To remove the artefact and the start empty part
                continue
            end
            % Adding corrections to meta feature
            if ftmark(nsbj,jj) == 0
                continue
            end
            ftnow = ft_sbj{nsbj}.ck{oridxft(jj)}.ft_ch(:,nft);    % The current feature values
            % ftnow (ftnow(:,1) == 0,:) = NaN;
            ftmat_this(nsbj,jj) = mean(ftnow,'omitnan');       % Average features across channels

        end

    end
    mean_allft(nft) = mean(ftmat_this(:),'omitnan');
    std_allft(nft) = std(ftmat_this(:),'omitnan');

end


%% Processing and save for FPCA analysis


ftall_mat_allnorm_strm = ftall_mat_allnorm(:,time_start_sampenough:end);
max_ck_real = size(ftall_mat_allnorm,2) - time_start_sampenough - epc_postasleep;

save('FPCADynamics.mat','ftall_mat_allnorm_strm','max_ck_real')

%% Save for other feature space analysis

% Clear up EEG to save space
clear eeg_asleepper_all
save("FeatureSpace_GroupAnalysis.mat")

%% Save global mean and std for normalisation

save ("GrpFt_MeanStd.mat",'mean_allft','std_allft')







