%% Supplementary Figure 9-11

% Author: JL

%% Supplementary Figure 9a

clear all
clc
load SuppFig9a.mat
onsetmark = epcs_tofit;
p_th = 0.025/length(pvec_dist);
%%%%%%%% Smoothed
figure
hold on
shadedErrorBar(tvec_now,xx_smooth,sdynam_ste(idxstartfit:end),'lineProps',{'-k','lineWidth',2})
plot(tvec_now,dd(1:end,1),'LineWidth',4,'Color','blue')
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('Sleep distance')
xlabel('Time (min)')
ax = gca;
ylim([0,7])
xlim([-30,10])
xticks([-30:10:10])
line([tvec_now(onsetmark),tvec_now(onsetmark)],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
% xlim([0,max_ck_real+epc_postasleep+1])
yl = ylim;
% 
% xticks([62:200:max_ck_real+epc_postasleep+1])
% xticklabels([-70:10:10])

for iii = 1:length(pvec_dist)
    tidx = iii+basetest_epc-1-idxstartfit;
    
    if pvec_dist(iii) < p_th   % Significant
        line([tvec_now(tidx),tvec_now(tidx+1)],[0.95*yl(2),0.95*yl(2)],'color','k','linewidth',1)
    end
end
scatter(tvec_now(idx_critic),dd(idx_critic,1),80,'red','filled')

%% Supplementary Figure 9b

% Warning: evaluation for bifurcation diagram takes very long

clear all
clc
load SuppFig9a.mat

c = dd(:,2);

param_min = params_optim;

r = param_min(1);
K = param_min(2);
h= param_min(3);
m = param_min(4);

x_theo = zeros(length(c));
c_3sol = [];

figure
hold on

syms x

for idxc = 1:length(c)
    
    % Roots of the theoretical equation when dx/dt = 0;
%     x = roots([1/10, -1, (c(idxc)+1/10), -1,0]);
    eqn = r*x*(1-x/K) - c(idxc)*x^2 / (x^2+h^2);
    xsol_all = solve(eqn==0);
    xsol_all = double(xsol_all);

    xsol = [];
    for jjj = 1:length(xsol_all)
        if isreal(xsol_all(jjj))
            if abs(xsol_all(jjj)) ~= 0
                xsol = [xsol,xsol_all(jjj)];
            end
        end
    end
    % Find 3 solution points
    if length(xsol) == 3
        c_3sol = [c_3sol,c(idxc)];
    end

    scatter(c(idxc)*ones(1,length(xsol)),xsol,2,'k')

end

xlabel('Control variable')
ylabel('System state')
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',3*get(gca,'ticklength'))
set(gca,'lineWidth',2)


%% Supplementary Figure 9c

clear all
clc
load SuppFig9c.mat

figure
edgesbin = [0:0.05:1];
pd = fitdist(r2all(mask_all),'Kernel');
pdfbin = 0:0.01:1;
ypdf = pdf(pd,pdfbin);

histogram(r2all(mask_all),edgesbin)
% histfit(tcrtc_new,30,'kernel')
hold on
plot(pdfbin,ypdf*(0.05*length(r2all(mask_all))),'LineWidth',2,'Color','r')
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('R2')
ylabel('Count')
xlim([0,1])
r2med = median(r2all(mask_all),'omitnan');
line([r2med,r2med],ylim,'LineStyle','--','LineWidth',2,'Color','r');

%% Supplementary Figure 10

clear all
clc
load SuppFig10.mat

% all subjects - SuppFig10a
sbjuniqs = unique(night_patidx);
numsbj = length(sbjuniqs);

mkrs = {"diamon","o","^","pentagram","v","square"}.';
maxlen1 = length(mkrs);
colorsall = distinguishable_colors(10);
colorsall = colorsall([2,1,4,5,9,10],:);

figure
hold on
clcidx = 1;
for ii = 1:numsbj

    idxthissbj = (night_patidx == sbjuniqs(ii));
    mkridx = floor(ii/maxlen1)+1;
    % if mkridx == 2 || mkridx==6 ||  mkridx == 3
        scatter(Y2all(idxthissbj,1),Y2all(idxthissbj,2),60,colorsall(clcidx,:),"filled",mkrs{mkridx})
    % else
    %     scatter(Y2all(idxthissbj,1),Y2all(idxthissbj,2),60,colorsall(clcidx,:),mkrs{mkridx})
    % end
    clcidx = clcidx+1;
    if clcidx>maxlen1
        clcidx = 1;
    end

end
legend
box off
set(gca, 'FontSize', 12);
set(gca, 'TickDir', 'out');
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

xlabel('TSNE-1')
ylabel('TSNE-2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consistent subject - SuppFig10b
sbjtest  =11;
night_thissbj = find(ismember(night_patidx,sbjtest));
Y2thissbj = Y2all(night_thissbj,:);

figure
scatter(Y2all(:,1),Y2all(:,2),'blue','filled')
hold on
scatter(Y2all(night_thissbj,1),Y2all(night_thissbj,2),'red','filled')
box off
set(gca, 'FontSize', 12);
set(gca, 'TickDir', 'out');
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

xlabel('TSNE-1')
ylabel('TSNE-2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inconsistent subject - SuppFig10c
sbjtest  = 6;
night_thissbj = find(ismember(night_patidx,sbjtest));
Y2thissbj = Y2all(night_thissbj,:);

figure
scatter(Y2all(:,1),Y2all(:,2),'blue','filled')
hold on
scatter(Y2all(night_thissbj,1),Y2all(night_thissbj,2),'red','filled')
box off
set(gca, 'FontSize', 12);
set(gca, 'TickDir', 'out');
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

xlabel('TSNE-1')
ylabel('TSNE-2')

%% Supplementary Figure 11a

clear all
clc
load SuppFig11a.mat

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
yline(criticalSValue, 'g--', 'LineWidth', 1.5)
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


%% Supplementary Figure 11b-c  

% Analysis codes for generating this result is inside the Process codes

clear all
clc
load SuppFig11bc.mat
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
ylabel('False alarm rate')
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)


figure
histogram(categorical(numFA_all(~idxmask)),'BarWidth',0.6)
box off
xlabel('False alarm')
ylabel('No. of nights')
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

mean(numFA_all(~idxmask),'omitnan')
std(numFA_all(~idxmask),'omitnan')






