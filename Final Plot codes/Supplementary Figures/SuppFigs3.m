%% Supplementary figure plots 3-4

% Author: JL

%% Supplementary Figure 3; 2D view of feature space

clear all
clc
load SuppFig3.mat

figure
histfit(dist_indiv_bedtime_Onset,20,'kernel')
box off
xlabel('Distance')
ylabel('No. of participants')
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlim([-1,20])

figure
gscatter(Y2d_tsne(:,1),Y2d_tsne(:,2),idxall,[],[],10)
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('TSNE-1')
ylabel('TSNE-2')










