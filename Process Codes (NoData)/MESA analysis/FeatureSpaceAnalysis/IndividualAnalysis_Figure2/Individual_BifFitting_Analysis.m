%% Analysis of the bifurcation function fittings

% Author: JL

%% Data loading

clear all
clc
load FitResults_Full.mat

fixed_tbl_results = tbl_results;

%% Correlation of R-squared to post sleep depth

num_sbj = size(fixed_tbl_results,1);
r2all = fixed_tbl_results.("R-squared");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sleep metrics
slp_scores = fixed_tbl_results.("Sleep-stages-trimmed");
post_slpecp = 20;

% Compute Post-sleep depth
depth_pos = NaN(num_sbj,1);

for ns = 1:num_sbj
    
    slpsc_now = slp_scores{ns};
    slpsc_pre_now = slpsc_now(1:end-post_slpecp);
    slpsc_post = slpsc_now(end-post_slpecp+1:end);

    depth_pos(ns) = mean(slpsc_post);

end

idxirr = (depth_pos>3);
depth_pos(idxirr) = [];
r2all(idxirr) = [];

% ifmonodec = logical(ones(num_sbj,1));
[r1,p1] = corr( depth_pos,r2all);

% Save for plots

save('ExtData_Fig1a.mat','depth_pos',"r2all")

%% Correlation between K distance and R2 Before exclusion

% Use K to approximate the attractor distance

num_sbj = size(fixed_tbl_results,1);

optim_params = fixed_tbl_results.("Optimal-Params");

K_all = NaN(num_sbj,1);
r2all = fixed_tbl_results.("R-squared");

for ns = 1:num_sbj
    
    K_all(ns) = optim_params{ns}(2);

end

ifexclude333 = zeros(num_sbj,1);
ifexclude333(K_all>15) = 1;    % Clean out one outlier

K_all_clean = K_all(~ifexclude333);
r2clean = r2all(~ifexclude333);

save('ExtData_Fig1b.mat','K_all_clean',"r2clean")


%% Clean up data, excludig participants with post-sleep consistency smaller than 1.5

num_sbj = size(fixed_tbl_results,1);
r2all = fixed_tbl_results.("R-squared");

slp_scores = fixed_tbl_results.("Sleep-stages-trimmed");
post_slpecp = 20;

ifpostdepthsmall = NaN(num_sbj,1);
postdepth_th = 1.5;
ifNoC = fixed_tbl_results.IfNoCritical;

for ns = 1:num_sbj
    
    slpsc_now = slp_scores{ns};
    slpsc_pre_now = slpsc_now(1:end-post_slpecp);
    slpsc_post = slpsc_now(end-post_slpecp+1:end);

    ifpostdepthsmall(ns) = (mean(slpsc_post)<postdepth_th);

end

sum(1-ifpostdepthsmall)
sum(ifNoC(logical(1-ifpostdepthsmall))==1)
sum(r2all(logical(1-ifpostdepthsmall))<0)

% Save for plots

ifexclude = ifpostdepthsmall;
ifexclude(ifNoC==1) = 1;
ifexclude(r2all<0) = 1;

r2all_new = r2all(logical(1-ifexclude));
mean(r2all_new)
std(r2all_new)
median(r2all_new)

save('Figure2h.mat',"r2all_new")
save('PostDepth.mat','ifpostdepthsmall')

%% Stats of tipping points & plot

t_crtic_all = fixed_tbl_results.("T-Critic");
t_crtic_new = t_crtic_all(logical(1-ifexclude));

mean(t_crtic_new)
std(t_crtic_new)
median(t_crtic_new)

figure
edgesbin = [-10:1:10];
pd = fitdist(t_crtic_new,'Kernel');
pdfbin = 0:0.01:1;
ypdf = pdf(pd,pdfbin);

histogram(t_crtic_new,edgesbin)
% histfit(tcrtc_new,30,'kernel')
hold on
% plot(pdfbin,ypdf*(0.05*length(t_crtic_new)),'LineWidth',2,'Color','r')
box off
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Tipping point time (min)')
ylabel('No. of participants')
% xlim([0,1])
r2med = median(t_crtic_new);
line([r2med,r2med],ylim,'LineStyle','--','LineWidth',2,'Color',[0.25,0.187,0.02]);
line([0,0],ylim,'LineStyle','--','LineWidth',2,'Color','r');

save('ExtData_Fig1c.mat','t_crtic_new')

%% Imports and evaluation of sleep latency and age

load epcsexc.mat    %Length of initial artefactual epochs to cut
load Sbjs_filtered_time.mat    % Raw PSG corresponding indexes
load("TS_SleepDistIndividuals.mat",'ifelig')    % The included subject IDs
load MESA_DemoInfo.mat    % Raw demo info
idx11 = epcs_exc_1(sbj_remained);
idx11 = idx11(logical(ifelig));

% Retrieve age info
ages_new = Age(sbj_remained);
ages_new = ages_new(logical(ifelig));

num_sbj = size(fixed_tbl_results,1);
r2all = fixed_tbl_results.("R-squared");
lat_all = NaN(num_sbj,1);
post_epc = 20;
slp_scores_ori = fixed_tbl_results.("Sleep-stages-original");

for ii = 1:num_sbj

    lat_all(ii) = (length(slp_scores_ori{ii})-20-idx11(ii)+1)/2;   % Latency in minutes
    
    if lat_all(ii) <=0
        lat_all(ii) = NaN;
    end
end

idxnan = ifexclude;
idxnan(isnan(lat_all)) = 1;

%% Take and test the relationship between sleep distance at bedtime (K) and R2 & Tipping points & latency after exclusion and cleaning

optim_params = fixed_tbl_results.("Optimal-Params");
t_crtic_all = fixed_tbl_results.("T-Critic");
r2all = fixed_tbl_results.("R-squared");

K_all = NaN(num_sbj,1);
slp_scores_ori = fixed_tbl_results.("Sleep-stages-original");

for ns = 1:num_sbj
    
    K_all(ns) = optim_params{ns}(2);

    % lat_all(ns) = (length(slp_scores_ori{ns})-20-idx11(ns))/2;   % Latency in minutes

end

ifexclude22 = idxnan;
ifexclude22(K_all>15) = 1;

K_all_clean = K_all(~ifexclude22);
t_crtic_clean = t_crtic_all(~ifexclude22);
lat_clean = lat_all(~ifexclude22);
r2clean = r2all(~ifexclude22);

save('ExtData_Fig1def.mat','K_all_clean','r2clean','lat_clean','t_crtic_clean')
save('Figure2i.mat','K_all_clean','t_crtic_clean')

[rr1,pp1] = corr(K_all_clean,(t_crtic_clean))
[rr2,pp2] = corr(K_all_clean,lat_clean)
[rr3,pp3] = corr((t_crtic_clean),lat_clean)
[rr4,pp4] = corr(r2clean,lat_clean)

%% Test at what sleep stage the tipping point happens

num_sbj = size(fixed_tbl_results,1);
r2all = fixed_tbl_results.("R-squared");

stage_tipp = NaN(num_sbj,1);

slp_scores = fixed_tbl_results.("Sleep-stages-trimmed");
t_crtic_all = fixed_tbl_results.("T-Critic");
tvec_all = fixed_tbl_results.("Time-vector");
hypno_exp = fixed_tbl_results.Expanded_Hyp;

t_crtic_all(logical(ifexclude)) = NaN;

for ii = 1:num_sbj

    if isnan(t_crtic_all(ii))   % Skip this
        continue
    end
    tidx = find(tvec_all{ii} == t_crtic_all(ii));   % Index of critical

    stage_tipp(ii) = hypno_exp{ii}(tidx);

end

% Quantify proportion
total_num = sum(1-isnan(stage_tipp))

disp('Proportion of Stage 0')
sum(stage_tipp==0)/total_num

disp('Proportion of Stage 1')
sum(stage_tipp==1)/total_num

disp('Proportion of Stage 2')
sum(stage_tipp>=2)/total_num

%% Test at what sleep stage the FPC2 peak happens  %% No exclusion at all, group stats

num_sbj = size(fixed_tbl_results,1);
r2all = fixed_tbl_results.("R-squared");

stage_tipp = NaN(num_sbj,1);

slp_scores = fixed_tbl_results.("Sleep-stages-trimmed");
t_crtic_all = fixed_tbl_results.("T-Critic");
tvec_all = fixed_tbl_results.("Time-vector");
hypno_exp = fixed_tbl_results.Expanded_Hyp;

for ii = 1:num_sbj

    if ifexclude(ii)  % Skip this
        continue
    end
    tidx = find(abs(tvec_all{ii}+ 4.50)<0.01);   % Index of critical

    if isempty(tidx)
        continue
    end

    stage_tipp(ii) = hypno_exp{ii}(tidx);

end

% Quantify proportion
total_num = sum(1-isnan(stage_tipp))

disp('Proportion of Stage 0')
sum(stage_tipp==0)/total_num

disp('Proportion of Stage 1')
sum(stage_tipp==1)/total_num

disp('Proportion of Stage 2')
sum(stage_tipp>=2)/total_num



%% Test correlations of fitting accuracy and tipping point to age % after exclusion

r2new = r2all(~ifexclude);
agesfilt = ages_new(~ifexclude);
t_crtic_all = fixed_tbl_results.("T-Critic");
tippfilt = t_crtic_all(~ifexclude);

idx_rem = logical(zeros(length(r2new),1));

[rr2,pp2] = corr(r2new(~idx_rem),agesfilt(~idx_rem))
[rr1,pp1] = corr(tippfilt,agesfilt)



