%% Extended Data Fig6

% JL

%% Extended Data Fig6a

clear all
clc
load ExtData_Fig6a.mat

figure('Position', [100, 100, 500, 400])
histogram(FPC1_var,20)
title('FPC1 Variance Explained Distribution')
xlabel('Variance Explained')
ylabel('Number of Subjects')
xlim([0,1])
ylim([0,150])
set(gca,'TickDir','out');
set(gca,'linewidth',2)
set(gca,'FontSize', 12)
set(gca,'ticklength',2*get(gca,'ticklength'))
box off


%% Extended Data Fig6b

clear all
clc
load ExtData_Fig6b.mat

figure
shadedErrorBar(1:size(bootstrap_mean_FPC1,2),bootstrap_mean_FPC1,abs(bootstrap_ste_FPC1),'lineProps',{'-k','lineWidth',2})
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',3*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('Magnitude (a.u)')
xlabel('Time (min)')
ax = gca;
xticks([0:200:size(bootstrap_mean_FPC1,2)+1])
xticklabels([-60:10:1])


%% Extended Data Fig6c

clear all
clc
load ExtData_Fig6c.mat

figure('Position', [100, 100, 500, 400])
bar(FPC1_prob)
title('Key FPC1 Probability Distribution')
xlabel('Feature Number')
ylabel('Probability')
set(gca,'TickDir','out');
set(gca,'linewidth',2)
set(gca,'FontSize', 12)
set(gca,'ticklength',2*get(gca,'ticklength'))
box off

%% Extended Data Fig6d

clear all
clc
load ExtData_Fig6d.mat

figure
edgesbin = [0:0.05:1];

histogram(rsq_thisft,edgesbin)
% histfit(tcrtc_new,30,'kernel')
hold on

box off
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Fitting accuracy (R2)')
ylabel('No. of participants')
xlim([0,1])
r2med = median(rsq_thisft,'omitmissing');
line([r2med,r2med],ylim,'LineStyle','--','LineWidth',2,'Color','r');


%% Extended Data Fig6e

clear all
clc
load ExtData_Fig6e.mat

figure
edgesbin = [0:0.05:1];
pd = fitdist(rsq_fpc_clean,'Kernel');
pdfbin = 0:0.01:1;
ypdf = pdf(pd,pdfbin);

histogram(rsq_fpc_clean,edgesbin)
% histfit(tcrtc_new,30,'kernel')
hold on
% plot(pdfbin,ypdf*(0.05*length(rsq_fpc_clean)),'LineWidth',2,'Color','r')
box off
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Fitting accuracy (R2)')
ylabel('No. of participants')
xlim([0,1])
r2med = median(rsq_fpc_clean);
line([r2med,r2med],ylim,'LineStyle','--','LineWidth',2,'Color','r');






