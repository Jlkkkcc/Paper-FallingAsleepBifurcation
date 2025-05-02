%% FPCA score analysis

% Junheng li

%%

clear all
clc
load FPCAscores.mat
load ftdscrp.mat
load FPCADynamics_Grp.mat

%% 

ft_use = [1:50];             % Feature number added;

ftall_mat_allnorm_strm = ftall_mat_allnorm_strm(ft_use,:);
ft_dscrp_add = ft_dscrp_add(ft_use);

%% Normalised FPCA scores

% Divide into positive and negative and separate them for discussion

%%%%%%%%%%%%%%%%%%%%%%%
% Positive FPC1 scores

pcs_posonly = pcscores(:,1);
idx_pos = pcs_posonly>0;
pcs_posonly(~idx_pos) = NaN;
pcs_posonly_norm = pcs_posonly./max(pcs_posonly);
ftidx_pos = find(pcs_posonly_norm>0.9)

%%%%%%%%%%%%%%%%%%%%%%%
% Negative FPC1 scores

pcs_negonly = pcscores(:,1);
idx_neg = pcs_negonly<0;
pcs_negonly(~idx_neg) = NaN;
pcs_negonly_norm = pcs_negonly./min(pcs_negonly);
ftidx_neg = find(pcs_negonly_norm>0.9)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patch back, normalisation to [-1,1]

pcnorm_all(idx_pos) = pcs_posonly_norm(idx_pos);
pcnorm_all(idx_neg) = -pcs_negonly_norm(idx_neg);
pcnorm_all = pcnorm_all';

%% Normalised FPCA scores for FPC2

% Divide into positive and negative and separate them for discussion

%%%%%%%%%%%%%%%%%%%%%%%
% Positive FPC2 scores

pcs_posonly = pcscores(:,2);
idx_pos = pcs_posonly>0;
pcs_posonly(~idx_pos) = NaN;
pcs_posonly_norm = pcs_posonly./max(pcs_posonly);
ftidx_pos = find(pcs_posonly_norm>0.9)

%%%%%%%%%%%%%%%%%%%%%%%
% Negative FPC2 scores

pcs_negonly = pcscores(:,2);
idx_neg = pcs_negonly<0;
pcs_negonly(~idx_neg) = NaN;
pcs_negonly_norm = pcs_negonly./min(pcs_negonly);
ftidx_neg = find(pcs_negonly_norm>0.9)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patch back, normalisation to [-1,1]
clear pcnorm_all
pcnorm_all(idx_pos) = pcs_posonly_norm(idx_pos);
pcnorm_all(idx_neg) = -pcs_negonly_norm(idx_neg);
pcnorm_all = pcnorm_all';










