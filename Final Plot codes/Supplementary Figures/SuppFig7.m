
%% Supplementary Figure 7
% JL


%% Supplementary Figure 7a

clear all
clc
load SuppFig7a.mat

stageNames = {'Awake', 'N1', 'N2', 'N3'};
f = figure;
f.Position(3:4) = [500,600];

tonset = find(tvec==0);

ax1=subplot(3,1,3);
plot(tvec,hyp_expand,'LineWidth',2)
ylim([0,3])
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
set(gca, 'YTick', 0:3, 'YTickLabel', stageNames);
line([0,0],ax1.YLim,'LineStyle','--','LineWidth',2,'Color','r');
hold on
crossingTimePoints = crossthisn.downwardCrossingTime;
crossingDetails = crossthisn;
if ~isempty(crossingTimePoints)
    for i = 1:length(crossingDetails)
        if crossingDetails(i).isRealCrossing

            tidx_cross = find(tvec>crossingDetails(i).downwardCrossingTime,1);
            h3 = plot(crossingDetails(i).downwardCrossingTime, hyp_expand(tidx_cross-1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Plot real crossing points
            % else
            %     h3 = plot(crossingDetails(i).downwardCrossingTime, criticalSValue, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Plot other crossing points
        end
    end
end

xlabel('Time (min)')
ax2=subplot(3,1,1:2);
plot(tvec,xx_smoothed,'k')
hold on
yline(criticalSValue, '--', 'LineWidth', 1.5,'Color',[0.929,0.694,0.125])
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('Sleep Distance')
% line([0,0],ax2.YLim,'LineStyle','--','LineWidth',2,'Color','r');

crossingTimePoints = crossthisn.downwardCrossingTime;
crossingDetails = crossthisn;
if ~isempty(crossingTimePoints)
    for i = 1:length(crossingDetails)
        if crossingDetails(i).isRealCrossing
            h3 = plot(crossingDetails(i).downwardCrossingTime, criticalSValue, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Plot real crossing points
            % else
            %     h3 = plot(crossingDetails(i).downwardCrossingTime, criticalSValue, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Plot other crossing points
        end
    end
end
line([0,0],ax2.YLim,'LineStyle','--','LineWidth',2,'Color','r');


%% Supplementary Figure 7b-c  

% Analysis codes for generating this result is inside the Process codes

clear all
clc
load SuppFig7bc.mat
figure
scatter(hypnofluct_alltrain,FA_rate_alltrain,'filled')
xlim([0.9,4])
ylim([-0.1,3.5])
p = polyfit(hypnofluct_alltrain,FA_rate_alltrain,1);
hold on
% xlim([-0.5,25])
xl = xlim;
dx = diff(xl)/100;
xvec = xl(1):dx:xl(2);
plot(xvec,xvec*p(1)+p(2),'r','LineWidth',2)

box off

xlabel('Fluctuation')
ylabel('Erroneous prediction rate')
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)


figure
histogram(categorical(numFA_all(~idxmask)),'BarWidth',0.6)
box off
xlabel('No. of erroneous prediction')
ylabel('No. of nights')
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

mean(numFA_all(~idxmask),'omitnan')
std(numFA_all(~idxmask),'omitnan')




