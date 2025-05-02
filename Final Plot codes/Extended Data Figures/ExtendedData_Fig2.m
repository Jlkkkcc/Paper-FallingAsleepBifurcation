%% Extended Data Figures plots Figure 2

% Author: Junheng Li

%% Extended Data Fig.2a-c

clear all
clc
load ExtData_Fig2abc.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frontal Fz-Cz

pvec_dist = pval_vec_Fz;
p_th = 0.025/length(pvec_dist);

onsetmark = epcs_tofit;
xxste = distall_ste1(time_start_sampenough:end);

figure
hold on
shadedErrorBar(tvec_now,xx_F,xxste(idxstartfit:end),'lineProps',{'-k','lineWidth',2})
plot(tvec_now,dd_F(1:end,1),'LineWidth',4,'Color','blue')
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('s')
xlabel('Time (min)')
ax = gca;
ylim([0,5])
xlim([-60,10])
xticks([-60:10:10])
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
% Add critical points
tcrtc_this = tcrtc_Fz;

scatter(tcrtc_this,dd_F((abs(tvec_now-tcrtc_this)<0.001),1),100,'red','filled','o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Posterior Cz-Oz

pvec_dist = pval_vec_Cz;
p_th = 0.025/length(pvec_dist);

onsetmark = epcs_tofit;
xxste = distall_ste2(time_start_sampenough:end);

figure
hold on
shadedErrorBar(tvec_now,xx_O,xxste(idxstartfit:end),'lineProps',{'-k','lineWidth',2})
plot(tvec_now,dd_O(1:end,1),'LineWidth',4,'Color','blue')
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('s')
xlabel('Time (min)')
ax = gca;
ylim([0,5])
xlim([-60,10])
xticks([-60:10:10])
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
% Add critical points
tcrtc_this = tcrtc_Cz;

scatter(tcrtc_this,dd_O((abs(tvec_now-tcrtc_this)<0.001),1),100,'red','filled','o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C4-M1

pvec_dist = pval_vec_C4;
p_th = 0.025/length(pvec_dist);

onsetmark = epcs_tofit;
xxste = distall_ste3(time_start_sampenough:end);

figure
hold on
shadedErrorBar(tvec_now,xx_C4,xxste(idxstartfit:end),'lineProps',{'-k','lineWidth',2})
plot(tvec_now,dd_C4(1:end,1),'LineWidth',4,'Color','blue')
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('s')
xlabel('Time (min)')
ax = gca;
ylim([0,5])
xlim([-60,10])
xticks([-60:10:10])
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
% Add critical points
tcrtc_this = tcrtc_C4;

scatter(tcrtc_this,dd_C4((abs(tvec_now-tcrtc_this)<0.001),1),100,'red','filled','o')


%% Extended Data Fig.2d-e

clear all
clc
load ExtData_Fig2de.mat

% Fz stats
r2all_new = r2all_Fz(logical(1-ifexclude));
mean(r2all_new)
std(r2all_new)
median(r2all_new)

figure
edgesbin = [0:0.05:1];
pd = fitdist(r2all_new,'Kernel');
pdfbin = 0:0.01:1;
ypdf = pdf(pd,pdfbin);
histogram(r2all_new,edgesbin)
% histfit(tcrtc_new,30,'kernel')
hold on
plot(pdfbin,ypdf*(0.05*length(r2all_new)),'LineWidth',2,'Color','r')
box off
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Fitting accuracy (R2)')
ylabel('No. of participants')
xlim([0,1])
r2med = median(r2all_new);
line([r2med,r2med],ylim,'LineStyle','--','LineWidth',2,'Color','r');

% Oz stats
clear r2all_new
r2all_new = r2all_Oz(logical(1-ifexclude));
mean(r2all_new)
std(r2all_new)
median(r2all_new)

figure
edgesbin = [0:0.05:1];
pd = fitdist(r2all_new,'Kernel');
pdfbin = 0:0.01:1;
ypdf = pdf(pd,pdfbin);
histogram(r2all_new,edgesbin)
% histfit(tcrtc_new,30,'kernel')
hold on
plot(pdfbin,ypdf*(0.05*length(r2all_new)),'LineWidth',2,'Color','r')
box off
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Fitting accuracy (R2)')
ylabel('No. of participants')
xlim([0,1])
r2med = median(r2all_new);
line([r2med,r2med],ylim,'LineStyle','--','LineWidth',2,'Color','r');




