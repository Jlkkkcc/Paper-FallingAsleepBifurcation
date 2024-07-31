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

% ifmonodec = logical(ones(num_sbj,1));
[r1,p1] = corr( depth_pos,r2all);

% Save for plots

save('Figure2G.mat','depth_pos',"r2all")

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

save('Figure2I.mat',"r2all_new")





