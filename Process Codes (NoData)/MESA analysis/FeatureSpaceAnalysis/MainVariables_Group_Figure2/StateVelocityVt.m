%% Feature space variable (group-level) -  state velocity v(t)

% Author: JL

%% Data loading

clear all
clc
load FeatureSpace_GroupAnalysis.mat
load GrpFt_MeanStd.mat

%% Computation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control parameters for how many features to use and how long time
% post-asleep to average on
ft_use = 1:47;             % Features used to evaluate distance
% ft_use = [10,11,15,41];
% ft_use = [13,15,45,41];
num_ft = length(ft_use);
N_samp_min = 200;

epc_postasleep = floor((10.5*60-epc_len)/jmp_len)+1;

max_epctotal = max(totalepc_sbj);

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
    % Feature matrix
    for jj = 1:max_epctotal

        if isnan(ftmark(nsbj,jj))            % To remove the artefact and the start empty part
            continue
        end
        if ftmark(nsbj,jj) % If not artefact
            ftnow = ft_sbj{nsbj}.ck{oridxft(jj)}.ft_ch(:,ft_use);    % The current feature values
            ftnow (ftnow(:,1) == 0,:) = NaN;
            ftmat_this(jj,:) = mean(ftnow,1,'omitnan');       % Average features across channels

            % Global normalisation
            ftmat_this_zsnorm(jj,:) = (ftmat_this(jj,:) - mean_allft) ./std_allft;
    
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Distance derivative
    
    for jj = 2:max_epctotal

        if isnan(ftmark(nsbj,jj))            % To remove the artefact and the start empty part
            continue
        end
        if isnan(ftmark(nsbj,jj-1))
            continue
        end

        ftvec_thistime = ftmat_this_zsnorm(jj-1,:);
        ftvec_next = ftmat_this_zsnorm(jj,:);
        if sum(isnan(ftvec_thistime))>0               % pdist currently cannot consider NaN
            continue
        end
        if sum(isnan(ftvec_next))>0               % pdist currently cannot consider NaN
            continue
        end
        % Euclidean distance
        dist_allsbjs_mat(nsbj,jj) = pdist2(ftvec_thistime,ftvec_next);

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrapping
% Previous lables cannot be used; So shuffling myself

numsamps_eachtime = sum((1-isnan(dist_allsbjs_mat)),1);
idx_enou_samp = find(numsamps_eachtime>=N_samp_min,1);
disp('Time to sleep onset in minutes:')
((max_epctotal-idx_enou_samp-epc_postasleep-1)*jmp_len+epc_len)/60

for jj = idx_enou_samp:max_epctotal

    num_samps_now = numsamps_eachtime(jj);
    distvecnow = dist_allsbjs_mat(:,jj);
    if num_samps_now<=N_samp_min
        idx_surr = find((1-isnan(distvecnow))== 1, num_samps_now);
    else
        idx_sbjall = find((1-isnan(distvecnow))== 1);
        ind_rand = randperm(length(idx_sbjall),N_samp_min);
        idx_surr = idx_sbjall(ind_rand);
    end
    
    dist_allsbjs_mat_bootstrap(idx_surr,jj) = dist_allsbjs_mat(idx_surr,jj);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check sample numbers
numsamp_tt = sum((1-isnan(dist_allsbjs_mat_bootstrap)),1);

distall_avg = mean(dist_allsbjs_mat_bootstrap,1,'omitnan');
distall_ste = std(dist_allsbjs_mat_bootstrap,[],1,'omitnan')./sqrt(numsamp_tt-1);

time_start_sampenough = find(numsamp_tt>=N_min,1);
max_ck_real = length(distall_avg) - time_start_sampenough - epc_postasleep;
tvec = -(max_ck_real-1)/20:0.05:10.5;

%% Stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stats
numsamp_tt = sum((1-isnan(distdiff_allsbjs_mat_boot)),1);

distall_avg_diff = mean(distdiff_allsbjs_mat_boot,1,'omitnan');
distall_ste_diff = std(distdiff_allsbjs_mat_boot,[],1,'omitnan')./sqrt(numsamp_tt-1);

time_start_sampenough = find(numsamp_tt>=N_samp_min,1);
max_ck_real = length(distall_avg_diff) - time_start_sampenough - epc_postasleep;
tvec = -(max_ck_real-1)/20:0.05:10.5;

base_pertest = 10*60;      % Window as baseline to test difference;
tail = 'both';       % Testing for decreasing

basetest_epc = floor((base_pertest - epc_len)/jmp_len)+1;

distmat_test = distdiff_allsbjs_mat_boot(:,time_start_sampenough:end);

pvec_dist = tvtest(distmat_test,basetest_epc,tail);

distall_avg_smooth = smoothdata(distall_avg,2,'movmedian',20);


%% Save for plots

save('Figure2B.mat','epc_len','jmp_len','max_ck_real','tvec','distall_avg_smooth','time_start_sampenough','distall_ste','pvec_dist','basetest_epc')








