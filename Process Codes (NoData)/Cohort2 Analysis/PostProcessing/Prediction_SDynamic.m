
%% Cosine similariy summary stats and plots

% Author: JL

%% Data loading

clear all
clc
% load SummaryTable.mat
load UpdatedSummaryTable.mat

%% Final plots (Leading to Figure4C)

nights_used = predictionSummaryTable.Number_of_Nights_Used;

colorsall = distinguishable_colors(10);
colorsall = colorsall([2,1,4,5,6,9,8],:);

cos_avg = predictionSummaryTable.AverageCosineSimilarity;

figure
violinplot(cos_avg,nights_used,'ViolinColor',colorsall)
box off
xlabel('No. of training nights')
ylabel('Cosine similarity')
set(gca, 'FontSize', 12);
set(gca, 'TickDir', 'out');
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

%% Statistics of Figure4C

idxunq = unique(nights_used);
lat_night= NaN(length(idxunq),2);
for ii = 1:length(idxunq)
    
    nightuse = idxunq(ii);
    nightidxnow = (nights_used==nightuse);

    latthisnight = cos_avg(nightidxnow);

    lat_night(ii,1) = mean(latthisnight,'omitmissing');
    lat_night(ii,2) = std(latthisnight,'omitmissing');
    

end

%% Estimation of increase in the cosine similarity score per additional training night

figure
errorbar(lat_night(:,1),lat_night(:,2),'-o','LineWidth',2,"MarkerSize",15,...
    "MarkerEdgeColor","blue","MarkerFaceColor","blue",'Color','k')
box off
xlabel('No. of training nights')
ylabel('Cosine similarity')
set(gca, 'FontSize', 12);
set(gca, 'TickDir', 'out');
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlim([0,8])
ylim([0.85,1.05])

p = polyfit([1:7],lat_night(:,1),1);
hold on
xl = xlim;
xplt = xl(1):0.01:xl(2);
plot(xplt,p(1)*xplt+p(2),'LineWidth',2,'Color','r','LineStyle','--')
xticks([1:7])








