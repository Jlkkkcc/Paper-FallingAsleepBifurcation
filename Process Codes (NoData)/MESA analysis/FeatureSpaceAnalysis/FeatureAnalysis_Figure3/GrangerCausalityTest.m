%% Causality test final for paper

% All methodology same as previous
% Author: JL

% Features used:
% FPCA1: Positive: 41, Negative: 15;    % For easy interpretation these two
% are changed but no effect in general
% FPCA2: Positive: 45, Negative: 13;

%% Data used
% This is different to other figures as this data is 6s sampled without
% overlapping;

clear all
clc

load FPCADynamics.mat

ftallmat11 = ftall_mat_allnorm_strm;
%% Testing for order of integration
% Normally the test should be KPSS and ADF test

epc_len = 6;
jmp_len = 3;
epc_postasleep = 209;

% ft_to_test = [13,33,45,47];
% ft_to_test = [13,15,6,41];
% ft_to_test = [13,15,45,41];
% ft_to_test = [6,13,41,43];
% ft_to_test = [6,13,46, 1];
ft_to_test = [41,46,6,17];

% ftall_mat_allnorm_strm = ftall_mat_allnorm(:,time_start_sampenough:end-epcs_to_asleep);
max_preepc = floor((60*60-epc_len)/jmp_len)+1;
ftall_mat_allnorm_strm = ftallmat11(:,end-epc_postasleep-max_preepc:end-epc_postasleep+19);

% ftall_mat_allnorm_strm = ftall_mat_allnorm(:,end-epc_postasleep-max_preepc:end);
% 
% [h,pValue,stat,cValue] = kpsstest(ftall_mat_allnorm_strm(ft_to_test,:))
% [h,pValue,stat,cValue] = adftest(ftall_mat_allnorm_strm(ft_to_test,:))

% Test the order of integration
% Nonstationary data can be repeatedly differenced into stationary thus to
% achieve the order of integration I(d)

Id_fts = zeros(length(ft_to_test),1);

for jj = 1:length(ft_to_test)

    ori_ts = ftall_mat_allnorm_strm(ft_to_test(jj),:);
    Id_fts(jj) = ord_int(ori_ts);

end

m = max(Id_fts);    % Maximum order of integration

disp(['Maxmimum integration: ',num2str(m)])

%% VAR model set-up: identifying optimal lag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selecting features to project to
ft_ids = ft_to_test;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lag_max = 30;
% Try to running on 4 top features only
h_LRtest = zeros(lag_max,1);
p_LRtest = zeros(lag_max,1);
num_params = zeros(lag_max,1);
aic_all = zeros(lag_max,1);          % Suitable models minimises AIC
bic_all = zeros(lag_max,1);

num_ft_var = 4;

for lag = 1:lag_max

    mdl_var = varm(num_ft_var,lag);

    [EstMdl,EstSE,logL,E] = estimate(mdl_var,ftall_mat_allnorm_strm(ft_ids,:)');
    results = summarize(EstMdl);
    numparams = results.NumEstimatedParameters;

    if lag == 1
        num_param_prev = numparams;
        LR_prev = logL;
        h_LRtest(lag) = NaN;
        p_LRtest(lag) = NaN;
    else
        [h_LRtest(lag),p_LRtest(lag)] = lratiotest(logL,LR_prev,numparams-num_param_prev);
        num_param_prev = numparams;
        LR_prev = logL;
    end
    aic_all(lag) = results.AIC;
    bic_all(lag) = results.BIC;

end

figure
plot(aic_all)
xlabel('Lags')
ylabel('AIC')

% Finding and building optimal lag VAR model based on AIC
lag_optm = find(aic_all == min(aic_all));

mdl_var = varm(num_ft_var,lag_optm);
[EstMdl_opt,EstSE,logL,E] = estimate(mdl_var,ftall_mat_allnorm_strm(ft_ids,:)');


%% Testing residual autocorrelations

[~,pValue1,~,~] = lbqtest(E(:,1));
[~,pValue2,~,~] = lbqtest(E(:,2));
[~,pValue3,~,~] = lbqtest(E(:,3));
[~,pValue4,~,~] = lbqtest(E(:,4));

p_th = 0.05/4;
h1 = pValue1<p_th;
h2 = pValue2<p_th;
h3 = pValue3<p_th;
h4 = pValue4<p_th;

% For each residual time-series test if there will be serial
% autocorrelation

% Test Durbin-Watson stats
dw_all = dwres(E');

if h1 || h2 || h3 || h4
    disp('Residual Serial Correlation detected using Ljung-Box test')
elseif max(dw_all-2) > 0.2
    disp('Residual Serial Correlation detected using DW test')
else
    disp('Can claim no significant serial correlation of residual')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Important! Increase model order p until the residual autocorrelation is
% not significantly correlated

% This step requires manual test; And in current status, optimal model
% order should +2 to result in no residual autocorrelation

% P-value used is 0.0125 for autocorrelation test

extra_lag = 2;
mdl_var = varm(num_ft_var,lag_optm+extra_lag);
[EstMdl_opt,EstSE,logL,E] = estimate(mdl_var,ftall_mat_allnorm_strm(ft_ids,:)');


%% Cointegration test

% Johansen cointegration test as MATLAB implemented

% Trace test & Max-eigenvalue test
[h_coi,pValue_coi,~,~] = jcitest(ftall_mat_allnorm_strm(ft_ids,:)',Test=["trace" "maxeig"]);

% Both results shows that C is full rank at 

% Testing using VEC model the cointegration rank
[~,C] = var2vec(EstMdl_opt.AR);
rank(C)

% If rank(C) == N_ft (Number of features) then it means the VAR model is
% stable at the level of yt.

%% GC test with extra coefficient

% Conduct the Chi-based Wald GC test on the first p-order: That is remove
% the last m order coefficients;

% GC test based on modified model
ts_totest = ftall_mat_allnorm_strm(ft_ids,:);
nft_totest = length(ft_ids);
vecids = 1:nft_totest;

adjmat_GCtest = zeros(nft_totest,nft_totest);
pmat_GCtest = zeros(nft_totest,nft_totest);
statmat_GCtest = zeros(nft_totest,nft_totest);
for ii = 1:nft_totest
    
    for jj = 1:nft_totest
        if jj == ii
            continue
        else
            vecids_now = vecids;
            vecids_now(vecids_now == ii) =[];
            vecids_now(vecids_now == jj) = [];
            % The VAR model should include the remaining two time-series as
            % a control for underlying model parameters
            [adjmat_GCtest(ii,jj),pmat_GCtest(ii,jj),statmat_GCtest(ii,jj)] = gctest(ts_totest(ii,:)',ts_totest(jj,:)',ts_totest(vecids_now,:)',Integration=m,NumLags=lag_optm+extra_lag);
%             [adjmat_GCtest(ii,jj),pmat_GCtest(ii,jj)] = gctest(ts_totest(ii,:)',ts_totest(jj,:)',Integration=m,NumLags=lag_optm);
        end
    end

end
ft_names = {ft_dscrp{ft_ids(:)}};

p_th = 0.05/2;

num_test = nft_totest*(nft_totest-1);
adjmat_GCtest_corBon = pmat_GCtest<(p_th);    % Bonferroni corrected
adjmat_GCtest_corBon = adjmat_GCtest_corBon - diag(diag(adjmat_GCtest_corBon));

G = digraph(adjmat_GCtest_corBon,ft_names);
p = plot(G,'MarkerSize',10,'LineWidth',1.5)

save('Figure3G.mat','p_th',"adjmat_GCtest_corBon",'pmat_GCtest','ft_names')

