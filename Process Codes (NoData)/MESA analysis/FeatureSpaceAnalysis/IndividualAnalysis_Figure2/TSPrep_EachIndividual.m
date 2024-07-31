%% Meta feature (Distance to sleep onset - EWS analysis)

% Author: JL

% This code is for the preparation of the sleep distance s(t) dynamics for
% each individual without artefact corrections

%% data loading

clear all
clc
load fteval_EWSallsbj_20Dec.mat
load GrpFt_MeanStd.mat

%% Parameters and pre-processing of state variable time-series

% Epoch definition
epc_len = 6;           % Epoch length in seconds
epc_lend = epc_len*Fs;                % Epoch length in data points
ovlap = 0.5;
jmp_len = epc_len - epc_len*ovlap;
jmp_lend = epc_lend - floor(epc_lend*ovlap);

% Artefact specific parameters
nrep = 1;
mul_rms = 2.5;  
plotart = 0;

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
    % Visualise artefact detection
end

%% Subject eligibility for EWS (artefacts checking)
% 04/02/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subject eligibility test
% Considering the Overlap, the length of artefact should be careful.

N_cont_art = 180;         % Maximum minutes length of continuous artefacts in second
max_artper = 0.5;           % Maximum percentage of artefacts
per_min = 10*60 / score_size;       % Minimum length of falling asleep period for EWS analysis
post_pernum = 21;     % Epochs of post-asleep (in scoring sizes)

time_ews = 30*60;    % Last N minutes for artefact filtering
time_ews_post = 10*60;   % Post asleep period for EWS calculation
jmp = epc_len*(1-ovlap);       % Time jump (considering overlap)

pre_gap = floor((time_ews-epc_len)/jmp)+1;
post_epc = floor((post_pernum*score_size-epc_len)/jmp)+1;
post_gap = floor((post_pernum*score_size-time_ews_post-epc_len)/jmp)+1;    % Number of epoch gap to the end of vector

ifelig = ones(1,num_sbjs);          % Marker of eligibility

for sbj = 1:num_sbjs
    
    % Time filtering
    if (length(oriscore_asleep{sbj})-post_pernum)<=per_min
        ifelig(sbj) = 0;
        continue
    end
    
    isart_sbj_now = isart_allsbj{sbj};
    % Artefact filtering only on last N minutes +1 minute post asleep
    len_asleep = (length(oriscore_asleep{sbj})-post_pernum)*score_size;
    if len_asleep >= time_ews
        isart_sbj_use = isart_sbj_now(end-pre_gap-post_epc+1:end-post_gap);
        totallen = time_ews+time_ews_post;
    else
        isart_sbj_use = isart_sbj_now(1:end-post_gap);
        totallen = len_asleep+time_ews_post;
    end

    [idx_st, num_cont] = artcont(isart_sbj_use);
    % Maximum continuous filtering
    max_cont = max(num_cont);
    max_contper = epc_len+jmp*(max_cont-1);      % Maximum continuous artefact length (s)
    if max_contper > N_cont_art
        ifelig(sbj) = 0;
        continue
    end
    % Maximum percentage filtering
    totalart_len = 0;
%     totallen = length(oriscore_asleep{sbj})*score_size;
    max_totalart = totallen*max_artper;
    for jj = 1:length(num_cont)
        
        totalart_len = totalart_len + epc_len+ (num_cont(jj)-1)*jmp;      % Total artefact length
    end
    if totalart_len >= max_totalart
        ifelig(sbj) = 0;
        continue
    end
    
end

%% Proc features

max_epctotal = max(totalepc_sbj);
epc_postasleep = floor((post_sleep_per*score_size-epc_len)/jmp)+1;        % Number of epochs for post-asleep period
ftmark = NaN(num_sbjs,max_epctotal);
ftidx_orid = NaN(num_sbjs,max_epctotal);

for sbj = 1:num_sbjs

    if isempty(eeg_asleepper_all{sbj})
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


%% Distance time-series construction revised

ft_list = [1:47];
% ft_list = [13,15,45,41];    % FPCA features

% Post-asleep time can be modified
% onset_per = 60;
% onsetepc = floor((onset_per-epc_len)/jmp)+1;
time_ews_post = 10*60;   % Post asleep period for EWS calculation

raw_post = 10*60;
post_gap = floor((post_pernum*score_size-raw_post-epc_len)/jmp)+1;    % Number of epoch gap to the end of vector
post_gap_ews = floor((post_pernum*score_size-time_ews_post-epc_len)/jmp)+1;    % Number of epoch gap to the end of vector

distmet = 'euclidean';    % Distance metric

ts_sbj_nomed = [];
% ft_no = [2,11,15,41];
max_len = (time_ews+time_ews_post-epc_len)/jmp+1;

% distmat_allsbj = NaN(num_sbjs,max_len);
distmat_allsbj_nomed = NaN(num_sbjs,max_len);

ns = 0;
num_allsbj = [];
% ftonset_all = [];

win_medfilt = 50;            % Window size for median filtering to replace NaN

for sbj = 1:num_sbjs
    if ifelig(sbj)
        ns = ns+1;
        isart_sbj_now = isart_allsbj{sbj};
        epclentotal = length(isart_sbj_now);
        len_asleep = (length(oriscore_asleep{sbj})-post_pernum)*score_size;
        if len_asleep >= time_ews
            isart_now = isart_sbj_now(end-pre_gap-post_epc+1:end-post_gap_ews);
            ck_idx = (epclentotal-pre_gap-post_epc+1:1:epclentotal-post_gap_ews);
            totallen = time_ews+time_ews_post;
%             onset_def = [time_ews:jmp:time_ews+onset_per-jmp]/jmp;
        else
            isart_now = isart_sbj_now(1:end-post_gap_ews);
            ck_idx = (1:1:epclentotal-post_gap_ews);
            totallen = len_asleep+time_ews_post;
%             onset_def = [len_asleep:jmp:len_asleep+onset_per-jmp]/jmp;
        end
        
        num_epcs_sbj = length(isart_now);
        num_allsbj(ns) = num_epcs_sbj;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sleep onset point evaluation
        
        onsetidx = (epclentotal - epc_postasleep-post_gap+1:epclentotal-post_gap);
        lengthonset = length(onsetidx);
        ftmatonset = NaN(lengthonset,length(ft_list));
        for ne = 1:lengthonset
            
            if isart_sbj_now(onsetidx(ne))
                ftmatonset(ne,:) = NaN(1,length(ft_list));
                continue
            end
            
            ftvec_allch = ft_sbj{sbj}.ck{onsetidx(ne)}.ft_ch(:,ft_list);
            ftvec_allch(ftvec_allch(:,1)==0,:) = NaN;
            ftmatonset(ne,:) = nanmean(ftvec_allch,1);
            % Global z-score normalisation
            ftmatonset(ne,:) = (ftmatonset(ne,:)-mean_allft(1:num_ft))./std_allft(1:num_ft);        % Global normlisation
            
        end
        ftonset_med = median(ftmatonset,1,'omitnan');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % No median filtering version
        ftmat = NaN(num_epcs_sbj,length(ft_list));
        
        for ne = 1:num_epcs_sbj
            
            ftvec_allch = ft_sbj{sbj}.ck{ck_idx(ne)}.ft_ch(:,ft_list);
            ftvec_allch(ftvec_allch(:,1)==0,:) = NaN;
            ftmat(ne,:) = mean(ftvec_allch,1,'omitnan');
            % Global z-score normalisation
            ftmat(ne,:) = (ftmat(ne,:)-mean_allft(1:num_ft))./std_allft(1:num_ft);        % Global normlisation
            
        end

        for ne = 1:num_epcs_sbj

            ts_sbj_nomed{sbj}.ftdist_noart(ne) = pdist2(ftmat(ne,:),ftonset_med,distmet);
            
        end
        distmat_allsbj_nomed(sbj,end-num_epcs_sbj+1:end) = ts_sbj_nomed{sbj}.ftdist_noart;
            
        stp_ftepcs_all_filt{sbj} = stp_ftepcs_all{sbj}(ck_idx(1:end));
    
    end

end
% 
clear eeg_asleepper_all ecg_now eeg_sbj_now
clear ft_sbj
save('TS_SleepDistIndividuals.mat')



