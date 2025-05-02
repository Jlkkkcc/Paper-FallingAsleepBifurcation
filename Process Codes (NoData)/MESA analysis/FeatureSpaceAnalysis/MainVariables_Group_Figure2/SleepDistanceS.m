
%% Group-level feature space variable computations - Sleep distance 

% Author: JL

%% Data loading
clear all
clc
load FeatureSpace_GroupAnalysis.mat
load GrpFt_MeanStd.mat

%% Computation of the Sleep distance s(t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control parameters for how many features to use and how long time
% post-asleep to average on
ft_use = 1:50;             % Features used to evaluate distance  % Before revision, 47 features

num_ft = length(ft_use);

epc_postasleep = floor((10.5*60-epc_len)/jmp_len)+1;
max_epctotal = max(totalepc_sbj);
% This controls how many epochs distance to compute on
time_todist = 10*60;         % Time of median point to compute
epcs_todist = floor((time_todist-epc_len)/jmp_len)+1;  
onset_mark = (max_epctotal - epc_postasleep+1):(max_epctotal - epc_postasleep+epcs_todist);

distall_avg = NaN(max_epctotal,1);
distall_ste = NaN(max_epctotal,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distance to sleep onset
% This is computed subject-wise

dist_allsbjs_mat = NaN(num_sbjs,max_epctotal);
dist_allsbjs_mat_bootstrap = NaN(num_sbjs,max_epctotal);

for nsbj = 1:num_sbjs
    
    ftmat_this = NaN(max_epctotal,num_ft);     % Empty feature matrix for this feature
    ftmat_this_zsnorm = NaN(max_epctotal,num_ft);
    oridxft = ftidx_orid(nsbj,:);     % The feature indexes of previous EWS calculation
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3-d feature matrix
    for jj = 1:max_epctotal

        if isnan(ftmark(nsbj,jj))            % To remove the artefact and the start empty part
            continue
        end
        if ftmark(nsbj,jj) == 0
            continue
        end
        ftnow = ft_sbj{nsbj}.ck{oridxft(jj)}.ft_ch(:,ft_use);    % The current feature values
        ftnow (ftnow(:,1) == 0,:) = NaN;
        ftmat_this(jj,:) = mean(ftnow,1,'omitnan');       % Average features across channels

        % Global normalisation
        ftmat_this_zsnorm(jj,:) = (ftmat_this(jj,:) - mean_allft) ./std_allft;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Distance to median of sleep onset period (marked)
    ftonset_this_norm = ftmat_this_zsnorm(onset_mark,:);
    ftonset_vec = median(ftonset_this_norm,1,'omitnan');
    
    for jj = 1:max_epctotal

        if isnan(ftmark(nsbj,jj))            % To remove the artefact and the start empty part
            continue
        end
        if ftmark(nsbj,jj) == 0
            continue
        end
        
        ftvec_thistime = ftmat_this_zsnorm(jj,:);
        if sum(isnan(ftvec_thistime))>0               % pdist currently cannot consider NaN
            continue
        end
        % Euclidean distance
        dist_allsbjs_mat(nsbj,jj) = pdist2(ftvec_thistime,ftonset_vec);

    end
    
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrapping
if size(dist_allsbjs_mat,2) ~= length(idx_bts)
    error('Check code')
end
for jj = 1:max_epctotal
    idx_bts_thistime = idx_bts{jj};
    if isempty(idx_bts_thistime)
        dist_allsbjs_mat_bootstrap(:,jj) = NaN;
    else
        dist_allsbjs_mat_bootstrap(idx_bts_thistime,jj) = dist_allsbjs_mat(idx_bts_thistime,jj);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check sample numbers
numsamp_tt = sum((1-isnan(dist_allsbjs_mat_bootstrap)),1);

distall_avg = mean(dist_allsbjs_mat_bootstrap,1,'omitnan');
distall_ste = std(dist_allsbjs_mat_bootstrap,[],1,'omitnan')./sqrt(numsamp_tt-1);

time_start_sampenough = find(numsamp_tt>=N_min,1);
max_ck_real = length(distall_avg) - time_start_sampenough - epc_postasleep;
tvec = -(max_ck_real-1)/20:0.05:10.5;

%% Run statistical testings and save

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stats
base_pertest = 10*60;      % Window as baseline to test difference;
tail = 'right';       % Testing for decreasing

basetest_epc = floor((base_pertest - epc_len)/jmp_len)+1;

distmat_test = dist_allsbjs_mat_bootstrap(:,time_start_sampenough:end);

pvec_dist = tvtest(distmat_test,basetest_epc,tail);

%% Save for bifurcation function fittings

clear ft_sbj
save('SleepDistance_Grp.mat')

%%

save('Figure2d.mat','distall_avg')






