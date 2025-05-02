
%% Individual subject FPCA analysis preparation

% JL

%% Data load

clear all
clc
load 'fteval_EWSallsbj_20Dec.mat'

%% Parameters and pre-processing of state variable time-series
% 03/02/2022
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
    if plotart
        plot(eeg_sbj_now.ch{ch})
        hold on
        for i = 1:num_epcs_max
            if isart_sbj(i)
                eeg_ch_epc = eeg_sbj_now.ch{ch}(stps_sbj(i):stps_sbj(i)+epc_lend-1);
                plot(stps_sbj(i):stps_sbj(i)+epc_lend-1,eeg_ch_epc,'r')
            end
        end
    end
end

%% Subject eligibility for model fitting (artefacts checking) - Last 60 minutes
% 04/02/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subject eligibility test
% Considering the Overlap, the length of artefact should be careful.

N_cont_art = 180;         % Maximum length of continuous artefacts in second
max_artper = 0.5;           % Maximum percentage of artefacts
per_min = 3*60 / score_size;       % Minimum length of falling asleep period for EWS analysis
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

%% Feature time series construction per subject (with no artefact corrections) - Same subjects as individual BIF fitting

% load GrpFt_MeanStd_Ftadd.mat

ft_list = [1:50];
% ft_list = [13,15,45,41];    % FPCA features
num_ft = length(ft_list);

% Post-asleep time can be modified
% onset_per = 60;
% onsetepc = floor((onset_per-epc_len)/jmp)+1;
time_ews_post = 10*60;   % Post asleep period for EWS calculation
time_ews = 60*60;    % Last N minutes for artefact filtering

pre_gap = floor((time_ews-epc_len)/jmp)+1;
post_gap = floor((post_pernum*score_size-time_ews_post-epc_len)/jmp)+1;    % Number of epoch gap to the end of vector

distmet = 'euclidean';    % Distance metric

ts_sbj_nomed = [];
% ft_no = [2,11,15,41];
max_len = (time_ews+time_ews_post-epc_len)/jmp+1;

% distmat_allsbj = NaN(num_sbjs,max_len);
distmat_allsbj_nomed = NaN(num_sbjs,max_len);

ns = 0;
num_allsbj = NaN(1,num_sbjs);
% ftonset_all = [];

win_medfilt = 50;            % Window size for median filtering to replace NaN

ft_fpca_all_mat = [];

for sbj = 1:num_sbjs

    if ~ifelig(sbj)
        continue
    end
    ns = ns+1;
    isart_sbj_now = isart_allsbj{sbj};
    epclentotal = length(isart_sbj_now);
    len_asleep = (length(oriscore_asleep{sbj})-post_pernum)*score_size;
    if len_asleep >= time_ews
        isart_now = isart_sbj_now(end-pre_gap-post_epc+1:end-post_gap);
        ck_idx = (epclentotal-pre_gap-post_epc+1:1:epclentotal-post_gap);
        totallen = time_ews+time_ews_post;
        %             onset_def = [time_ews:jmp:time_ews+onset_per-jmp]/jmp;
    else
        isart_now = isart_sbj_now(1:end-post_gap);
        ck_idx = (1:1:epclentotal-post_gap);
        totallen = len_asleep+time_ews_post;
        %             onset_def = [len_asleep:jmp:len_asleep+onset_per-jmp]/jmp;
    end

    num_epcs_sbj = length(isart_now);
    num_allsbj(sbj) = num_epcs_sbj;

    ftmat_thissbj = NaN(num_ft,num_epcs_sbj);


    for nft = 1:length(ft_list)
        
        ftid =  ft_list(nft);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sleep onset point evaluation

        onsetidx = (epclentotal - epc_postasleep-post_gap+1:epclentotal-post_gap);
        lengthonset = length(onsetidx);
        ftmatonset = NaN(lengthonset,1);
        for ne = 1:lengthonset
            thisartidx = find(ck_idx==onsetidx(ne));

            if isart_now(thisartidx)
                ftmatonset(ne,:) = NaN;
                continue
            end

            ftvec_allch = ft_sbj_patch{sbj}.ck{onsetidx(ne)}.ft_ch(:,ftid);
            % ftvec_allch(ftvec_allch(:,1)==0,:) = NaN;
            ftmatonset(ne,:) = nanmean(ftvec_allch,1);
            % Global z-score normalisation
            ftmatonset(ne,:) = (ftmatonset(ne,:)-mean_allft(ftid))./std_allft(ftid);        % Global normlisation

        end
        ftonset_med = median(ftmatonset,1,'omitnan');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % No median filtering version
        ftmat = NaN(num_epcs_sbj,1);

        for ne = 1:num_epcs_sbj

            ftvec_allch = ft_sbj_patch{sbj}.ck{ck_idx(ne)}.ft_ch(:,ftid);
            % ftvec_allch(ftvec_allch(:,1)==0,:) = NaN;
            ftmat(ne,:) = mean(ftvec_allch,1,'omitnan');
            % Global z-score normalisation
            ftmat(ne,:) = (ftmat(ne,:)-mean_allft(ftid))./std_allft(ftid);        % Global normlisation

        end
        %
        tsnow = NaN(1,num_epcs_sbj);
        for ne = 1:num_epcs_sbj

            tsnow(ne) = pdist2(ftmat(ne,:),ftonset_med,distmet);

        end

        ftmat_thissbj(nft,:) = tsnow;


       

    end
    ft_fpca_all_mat{sbj} = ftmat_thissbj;

end


%% Processing

ft_dscrp = [ft_dscrp ft_dscrp_new(5:7)];
ftall_mat_allnorm_strm_all_subj = {};
idx_subj_used = [];
max_ck_real_all_subj = [];
for i = 1:length(ft_fpca_all_mat)
    if ~isempty(ft_fpca_all_mat{1,i})
        
        ftall_mat_allnorm_strm_subj = ft_fpca_all_mat{1,i};
        max_ck_real_subj = length(ft_fpca_all_mat{1,i})-200;

        % 1) Find indices of columns containing any NaN
        badColMask = any(isnan(ftall_mat_allnorm_strm_subj), 1);
        badCols    = find(badColMask);

        % 2) Adjust max_ck_real_subj if bad column index is larger
        %    than max_ck_real_subj
        for cIdx = 1 : length(badCols)
            if badCols(cIdx) > max_ck_real_subj
                max_ck_real_subj = max_ck_real_subj + 1;
            end
        end

        % 3) Remove the bad columns
        ftall_mat_allnorm_strm_subj(:, badCols) = [];

        max_ck_real_all_subj = [max_ck_real_all_subj max_ck_real_subj];
        ftall_mat_allnorm_strm_all_subj{end+1} = ftall_mat_allnorm_strm_subj;
        idx_subj_used = [idx_subj_used i];
    end
end


clear eeg_asleepper_all
save('Ft50Dynamics_IndivFPCA.mat', 'max_ck_real_all_subj', ...
    'ftall_mat_allnorm_strm_all_subj','idx_subj_used')

save('ft_dscrp.mat','ft_dscrp');
























