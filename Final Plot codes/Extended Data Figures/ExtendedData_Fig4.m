%% Extended Data Fig4

% JL

%% Extended Data Fig4a-i

clear all
clc
load ExtData_Fig4atoi.mat

% The features plotted in Extended Data Fig4a-i; To plot specific features,
% please change the plotting index "ft_to_plot" to relevant features
% indices; Refer to description (ft_dscrp) to match feature indexes with names
% Notice that, the y-axis scales need to be changed for different features
% to adjust to their ranges.

% Features plotted in the figure:
% Theta band power: 2; theta-beta ratio: 9; Dwelling time: 33; Spectral
% Slope: 50; Spectral centroid: 46; Prediction Error: 47; LZ complexity:
% 49; Peak theta frequency:13; Alpha temporal coherence: 18;

ft_to_plot = 2;

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

% ylim([-0.4,1.6])
% yticks([-0.4:0.4:1.6])
line([0,0],ylim,'Color','r','LineStyle','--','LineWidth',1.5)
xlim([-60,10])
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
title(ft_dscrp_add{ft_to_plot})

%% Extended Data Fig4j

clear all
clc
load ExtData_Fig4j.mat

figure
hold on
plot(tvec,xx_smooth/4,'Color','k','LineWidth',3)
plot(tvec,dd(1:end,1)/4,'LineWidth',3,'Color','blue')
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('Magnitude (a.u.)')
xlabel('Time (min)')
ax = gca;
ylim([0,1])
xlim([-60,0])
xticks([-60:10:0])

yl = ylim;

% Tipping point
scatter(t_critc,dd(find(tvec==t_critc),1)/4,80,'filled','MarkerFaceColor','r')

%% Extended Data Fig4k

clear all
clc
load ExtData_Fig4k.mat

figure
hold on

scatter(x2,x1,80,tvec,'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Magnitude FPC2 (a.u.)')
ylabel('Magnitude FPC1 (a.u.)')
c = colorbar;
c.Limits = [-60,0];
c.Ticks = [-60:10:0];
c.Label.String = 'Time (min)';


%% Extended Data Fig4l

clear all
clc
load ExtData_Fig4l.mat


% ftall_mat_allnorm_strm = ftall_mat_allnorm(:,time_start_sampenough:end);
tvec = -(max_ck_real-1)/20:0.05:10.5;
% nt_real = size(ftall_mat_allnorm_strm,2);
% slonset_n = nt_real - epc_postasleep-1;
% tvec = (-slonset_n*epc_len)/60:epc_len/60:(epc_postasleep*epc_len)/60;

figure
plot(tvec, ftall_mat_allnorm_strm(1,:))
xlabel('Time (min)')
ylabel('z-score')
hold on
line(xlim,[0 0],'Color','k','LineStyle','--','LineWidth',1.5)

ylim([14,22])
% yticks([-0.4:0.4:1.6])
line([0,0],ylim,'Color','r','LineStyle','--','LineWidth',1.5)
xlim([-60,10])
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)


