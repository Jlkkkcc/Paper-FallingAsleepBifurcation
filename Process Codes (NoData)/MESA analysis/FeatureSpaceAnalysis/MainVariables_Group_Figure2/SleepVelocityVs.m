
%% Evaluation Group-level sleep velocity v_s 

% Author: JL

%% Data loading

clear all
clc
load FeatureSpace_GroupAnalysis.mat
load GrpFt_MeanStd.mat

%% Evaluation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control parameters for how many features to use and how long time
% post-asleep to average on
ft_use = 1:47;             % Features used to evaluate distance
% ft_use = [10,11,15,41];
% ft_use = [13,15,45,41];
num_ft = length(ft_use);

max_epctotal = max(totalepc_sbj);
% This controls how many epochs distance to compute on
time_todist = 10*60;         % Time of median point to compute
epcs_todist = floor((time_todist-epc_len)/jmp_len)+1;
epc_postasleep = floor((10.5*60-epc_len)/jmp_len)+1;
onset_mark = (max_epctotal - epc_postasleep+1):(max_epctotal - epc_postasleep+epcs_todist);

distall_avg = NaN(max_epctotal,1);
distall_ste = NaN(max_epctotal,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distance to sleep onset
% This is computed subject-wise

theta_allsbjs_mat = NaN(num_sbjs,max_epctotal);
theta_allsbjs_mat_bootstrap = NaN(num_sbjs,max_epctotal);

dotp_allsbjs_mat = NaN(num_sbjs,max_epctotal);
dotp_allsbjs_mat_bootstrap = NaN(num_sbjs,max_epctotal);

runwin = 20;

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
        % Adding corrections to meta feature
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Phase angles between current direction (from t to t+1) and
    % direction towards sleep onset cluster

    for jj = 1:max_epctotal-runwin-1

        if isnan(ftmark(nsbj,jj))            % To remove the artefact and the start empty part
            continue
        end
        if isnan(ftmark(nsbj,jj+1))            % To remove the artefact and the start empty part
            continue
        end
        % Adding corrections to meta feature
        if ftmark(nsbj,jj) == 0
            continue
        end
        % Adding corrections to meta feature
        if ftmark(nsbj,jj+1) == 0
            continue
        end
        
        ftvec_thistime = ftmat_this_zsnorm(jj:(jj+runwin),:);
        ftvec_next = ftmat_this_zsnorm((jj+1):(jj+1+runwin),:);
        % if sum(isnan(ftvec_thistime))>0               % pdist currently cannot consider NaN
        %     continue
        % end
        % if sum(isnan(ftvec_thistime))>0               % pdist currently cannot consider NaN
        %     continue
        % end
        ftvec_thistime = mean(ftvec_thistime,1,'omitmissing');
        ftvec_next = mean(ftvec_next,1,'omitmissing');
        
        vecnext = ftvec_next - ftvec_thistime;
        veconset = ftonset_vec - ftvec_thistime;
        veconset_next = ftonset_vec - ftvec_next ;

        % Calculate phase angles
        vecnext_norm = vecnext./sqrt(vecnext*vecnext');
        veconset_norm = veconset./sqrt(veconset*veconset');
        theta = acos(vecnext_norm*veconset_norm');

        d1 = norm(veconset);
        d2 = norm(veconset_next);
        dist_allsbjs_mat(nsbj,jj) = (d1-d2);

        theta_allsbjs_mat(nsbj,jj) = theta;
        dotp_allsbjs_mat(nsbj,jj) = sqrt(vecnext*vecnext') * cos(theta);    % Projections onto Onset vector
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrapping
% Previous lables cannot be used; So shuffling myself

max_epctotal = size(amp_allsbjs_mat,2) - runwin;
numsamps_eachtime = sum((1-isnan(dotp_allsbjs_mat)),1);
idx_enou_samp = find(numsamps_eachtime>=N_min,1);
disp('Time to sleep onset in minutes:')
((max_epctotal-idx_enou_samp-epc_postasleep-1)*jmp_len+epc_len)/60

dotp_avg = NaN(1,max_epctotal);
dotp_ste = NaN(1,max_epctotal);

for jj = idx_enou_samp:max_epctotal-1

    num_samps_now = numsamps_eachtime(jj);
    distvecnow = theta_allsbjs_mat(:,jj);
    if num_samps_now<=N_min
        idx_surr = find((1-isnan(distvecnow))== 1, num_samps_now);
    else
        idx_sbjall = find((1-isnan(distvecnow))== 1);
        ind_rand = randperm(length(idx_sbjall),N_min);
        idx_surr = idx_sbjall(ind_rand);
    end
       
    dotp_allsbjs_mat_bootstrap(idx_surr,jj) = dotp_allsbjs_mat(idx_surr,jj);

    dotp_now = dotp_allsbjs_mat(idx_surr,jj);
    dotp_avg(jj) = mean(dotp_now);
    dotp_ste(jj) = std(dotp_now)./sqrt(num_samps_now-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check sample numbers
numsamp_tt = sum((1-isnan(theta_allsbjs_mat_bootstrap)),1);

% distall_avg = mean(dist_allsbjs_mat_bootstrap,1,'omitnan');
% distall_ste = std(dist_allsbjs_mat_bootstrap,[],1,'omitnan')./sqrt(numsamp_tt-1);
time_start_sampenough = find(numsamp_tt>=N_min,1);
max_ck_real = max_epctotal - time_start_sampenough - epc_postasleep;
tvec = -(max_ck_real-1)/20:0.05:10.5;

%% Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stats
base_pertest = 10*60;      % Window as baseline to test difference;
tail = 'left';       % Testing for decreasing

basetest_epc = floor((base_pertest - epc_len)/jmp_len)+1;

distmat_test = dotp_allsbjs_mat_bootstrap(:,time_start_sampenough:end);

pvec_dist = tvtest(distmat_test,basetest_epc,tail);

tvec = -(max_ck_real-1)/20:0.05:10.5;
p_th = 0.025 /(length(pvec_dist)-20);

%% Save for plots

save('Figure2C.mat','epc_len','jmp_len','max_ck_real','tvec','dotp_avg','dotp_ste','time_start_sampenough','pvec_dist','basetest_epc')
































