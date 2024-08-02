

%% Plot individual features (Figure 3B and 3E)

clear all
clc
% load ftdynam_6ssampOvlap3s_Onset1minCorrAlign_N200_FPCA_revMay2024.mat
load FPCADynamics.mat
load ftdscrp.mat   % The feature code names to match the index

%% Individual feature time-series plots-46

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change this to match the feature you want to plot
% The four key feature indexes used in the paper:
% FPC1: Total EEG power (41), Spectrum centroid (46)
% FPC2: Delta-to-Alpha ratio (6), Theta temporal coherence (17)

% Refer to description (ft_dscrp) to match feature indexes with names
ft_to_plot = 46;

% ftall_mat_allnorm_strm = ftall_mat_allnorm(:,time_start_sampenough:end);
tvec = -(max_ck_real-1)/20:0.05:10.5;
% nt_real = size(ftall_mat_allnorm_strm,2);
% slonset_n = nt_real - epc_postasleep-1;
% tvec = (-slonset_n*epc_len)/60:epc_len/60:(epc_postasleep*epc_len)/60;

figure
plot(tvec, ftall_mat_allnorm_strm(ft_to_plot,:))
xlabel('Time (min)')
ylabel('z-score')
hold on
line(xlim,[0 0],'Color','k','LineStyle','--','LineWidth',1.5)
line([0,0],ylim,'Color','r','LineStyle','--','LineWidth',1.5)

ylim([-0.5,2])
yticks([-0.5:0.5:2])
xlim([-60,10])
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

% Check y-lims

%% Individual feature time-series plots -41

ft_to_plot = 41;

% ftall_mat_allnorm_strm = ftall_mat_allnorm(:,time_start_sampenough:end);
tvec = -(max_ck_real-1)/20:0.05:10.5;
% nt_real = size(ftall_mat_allnorm_strm,2);
% slonset_n = nt_real - epc_postasleep-1;
% tvec = (-slonset_n*epc_len)/60:epc_len/60:(epc_postasleep*epc_len)/60;

figure
plot(tvec, ftall_mat_allnorm_strm(ft_to_plot,:))
xlabel('Time (min)')
ylabel('z-score')
hold on
line(xlim,[0 0],'Color','k','LineStyle','--','LineWidth',1.5)
line([0,0],ylim,'Color','r','LineStyle','--','LineWidth',1.5)

ylim([-2,0.5])
yticks([-2:0.5:0.5])
xlim([-60,10])
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

% Check y-lims

%% Individual feature time-series plots - 6

ft_to_plot = 6;

% ftall_mat_allnorm_strm = ftall_mat_allnorm(:,time_start_sampenough:end);
tvec = -(max_ck_real-1)/20:0.05:10.5;
% nt_real = size(ftall_mat_allnorm_strm,2);
% slonset_n = nt_real - epc_postasleep-1;
% tvec = (-slonset_n*epc_len)/60:epc_len/60:(epc_postasleep*epc_len)/60;

figure
plot(tvec, ftall_mat_allnorm_strm(ft_to_plot,:))
xlabel('Time (min)')
ylabel('z-score')
hold on
line(xlim,[0 0],'Color','k','LineStyle','--','LineWidth',1.5)
line([0,0],ylim,'Color','r','LineStyle','--','LineWidth',1.5)

ylim([-1,0.4])
yticks([-1:0.2:0.4])
xlim([-60,10])
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

% Check y-lims
%% Individual feature time-series plots - 17

ft_to_plot = 17;

% ftall_mat_allnorm_strm = ftall_mat_allnorm(:,time_start_sampenough:end);
tvec = -(max_ck_real-1)/20:0.05:10.5;
% nt_real = size(ftall_mat_allnorm_strm,2);
% slonset_n = nt_real - epc_postasleep-1;
% tvec = (-slonset_n*epc_len)/60:epc_len/60:(epc_postasleep*epc_len)/60;

figure
plot(tvec, ftall_mat_allnorm_strm(ft_to_plot,:))
xlabel('Time (min)')
ylabel('z-score')
hold on
line(xlim,[0 0],'Color','k','LineStyle','--','LineWidth',1.5)
line([0,0],ylim,'Color','r','LineStyle','--','LineWidth',1.5)

ylim([-0.2,1])
yticks([-0.2:0.2:1])
xlim([-60,10])
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

% Check y-lims
