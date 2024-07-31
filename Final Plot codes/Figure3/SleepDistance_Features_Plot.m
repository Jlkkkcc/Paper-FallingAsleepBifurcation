%% Bifurcation diagram - Sleep distance s(t) versus the feature dynamics

% Author: JL

%% Load data

clear all
clc

load Figure2A.mat
load FPCADynamics.mat

%% Bifurcation diagram against all types of control parameters as features

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change this to match the feature you want to plot
% The four key feature indexes used in the paper:
% FPC1: Total EEG power (41), Spectrum centroid (46)
% FPC2: Delta-to-Alpha ratio (6), Theta temporal coherence (17)

% Refer to description (ft_dscrp) to match feature indexes with names for
% other features
ctr_ft = 17;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tvec = -(max_ck_real-1)/20:0.05:10.5;
tvec = tvec(idxstartfit:end);

% tvec = -(max_ck_real-2)/20:0.05:10.5;
ftplt_smooth = smoothdata(ftall_mat_allnorm_strm(ctr_ft,idxstartfit:end),2,'movmedian',10);
xxsmthplt = smoothdata(xx,2,'movmedian',10);

figure 
hold on
scatter(ftplt_smooth,xxsmthplt,20,tvec,'filled','Marker','o')
% scatter(ftall_mat_allnorm_strm(ctr_ft,:),xx)

colorbar

xlabel('z-score')
ylabel('Sleep Distance')
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
% yticks()

zero_dix = find(tvec == 0);
pltrange = zero_dix:zero_dix+20;
scatter(ftplt_smooth(pltrange),xxsmthplt(pltrange),20,'r','filled')

% Modify colorbar
gg = parula(length(tvec));
gg(pltrange,:) = repmat([1,0,0],length(pltrange),1);
colormap(gg)
ylim([0,5])
yticks([0:5])

% clabel('Time (min)')

% figure 
% plot(smoothdata(ftall_mat_allnorm_strm(ctr_ft,idxstartfit:end),2,'movmedian',10),smoothdata(xx,2,'movmedian',10),'-k')
% % scatter(ftall_mat_allnorm_strm(ctr_ft,:),xx)

%% Testing linear relationships between s and features

[r,p] = corr(xxsmthplt',ftplt_smooth')





