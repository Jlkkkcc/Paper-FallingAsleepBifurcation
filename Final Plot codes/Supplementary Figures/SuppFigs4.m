%% Supplementary Figure 4: Confounding test for critical slowing down dynamics

clear all
clc
load SuppFig4.mat

p_th = 0.025 / length(pvec_slp);         % Bonferroni correction
% p_th = 0.025;

xl = [-50,-29];
xticksthis = [-50:5:-30];

idxplt = find(tews_plt<-29);
tews_plt = tews_plt(idxplt);

% onsetmark = find(tews_plt == 0);      % Sleep onset marker

lengap = length(tews) - length(tews_plt);

f = figure;
f.Position(3:4) = [550,650];
hold on

avg_score1 = mean(avg_score,1,'omitnan');
std_score1 = std(avg_score,[],1,'omitnan')./sqrt(kk-1);

avg_score1 = avg_score1(end-length(tews)+1:end);
std_score1 = std_score1(end-length(tews)+1:end);

avg_score1 = avg_score1(idxplt);
std_score1 = std_score1(idxplt);

subplot(3,1,1)
shadedErrorBar(tews_plt, avg_score1,std_score1,'lineProps',{'k','lineWidth',2})
title('Avg. sleep stages')
ylabel('Sleep stages')
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylim([-0.5,1])
yl = ylim;
% line([tews_plt(onsetmark),tews_plt(onsetmark)],yl,'LineStyle','--','LineWidth',2,'Color','r')
xlim(xl)
xticks(xticksthis)
for iii = 1:length(pvec_slp)
    tidx = iii+base_epcs-1-lengap;
    
    if pvec_slp(iii) < p_th   % Significant
        line([tews_plt(tidx),tews_plt(tidx+1)],[2.3,2.3],'color','k','linewidth',1)
    end
end

allewsts =  ews_all{1};
ewsavg = median(allewsts,1,'omitnan');
% ewsavg = [NaN(1,length(t)-length(ewsavg)),ewsavg];
stdews = std(allewsts,[],1,'omitnan')./sqrt(kk-1);
% stdews = [NaN(1,length(t)-length(stdews)),stdews];

ewsavg = ewsavg(idxplt);
stdews = stdews(idxplt);

subplot(3,1,2)
shadedErrorBar(tews_plt, ewsavg,stdews,'lineProps',{'k','lineWidth',2})
% title(ews_to_avg{1})
xlim(xl)
xticks(xticksthis)
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylim([0.35,0.45])
yl = ylim;
ylabel('Autocorrelation')
% line([tews_plt(onsetmark),tews_plt(onsetmark)],yl,'LineStyle','--','LineWidth',2,'Color','r')
pvecnow = pvec_ews{1};
for iii = 1:length(pvecnow)
    tidx = iii+base_epcs-1-lengap;
    
    if pvecnow(iii) < p_th   % Significant
        line([tews_plt(tidx),tews_plt(tidx+1)],[0.99*yl(2),0.99*yl(2)],'color','k','linewidth',1)
    end
end


allewsts =  ews_all{2};
ewsavg = median(allewsts,1,'omitnan');
% ewsavg = [NaN(1,length(t)-length(ewsavg)),ewsavg];

stdews = std(allewsts,[],1,'omitnan')./sqrt(kk-1);
% stdews = [NaN(1,length(t)-length(stdews)),stdews];

ewsavg = ewsavg(idxplt);
stdews = stdews(idxplt);

subplot(3,1,3)
shadedErrorBar(tews_plt, ewsavg,stdews,'lineProps',{'k','lineWidth',2})
% title(ews_to_avg{2})
xlim(xl)
xticks(xticksthis)
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('Variance')
ylim([1.1,1.5])
yl = ylim;
% line([tews_plt(onsetmark),tews_plt(onsetmark)],yl,'LineStyle','--','LineWidth',2,'Color','r')
pvecnow = pvec_ews{2};
for iii = 1:length(pvecnow)
    tidx = iii+base_epcs-1-lengap;
    
    if pvecnow(iii) < p_th   % Significant
        line([tews_plt(tidx),tews_plt(tidx+1)],[0.99*yl(2),0.99*yl(2)],'color','k','linewidth',1)
    end
end

xlabel('Time (min)')


