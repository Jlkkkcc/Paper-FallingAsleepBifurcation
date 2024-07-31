%% All final plots for Figure 2

% Author: JL

% Data were processed and analysed in MESA analysis folder; Only final data
% are shown and included here for plotting with no process data.

%% Figure 2A

clear all
clc
load Figure2A.mat

p_th = 0.025 / length(pvec_dist);
onsetmark = epcs_tofit;
xxste = distall_ste(time_start_sampenough:end);
%%%%%%%% Smoothed
figure
hold on
shadedErrorBar(tvec_now,xx_smooth,xxste(idxstartfit:end),'lineProps',{'-k','lineWidth',2})
plot(tvec_now,dd(1:end,1),'LineWidth',4,'Color','blue')
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
% Tipping point
scatter(t_critc,dd(find(tvec_now==t_critc),1),80,'filled','MarkerFaceColor','r')


%% Figure 2B

clear all
clc
load Figure2B.mat

timeplot = 60*60;
epcs_toplot = floor((timeplot-epc_len)/jmp_len)+1;

idxstartplot = max_ck_real-epcs_toplot;

tvec_new = tvec(idxstartplot:end);
distvec_smooth = distall_avg_smooth(time_start_sampenough:end);
distvec_smooth_new = distvec_smooth(idxstartplot:end);

distste = distall_ste(time_start_sampenough:end);
distste_new = distste(idxstartplot:end);

onset_new = max_ck_real - idxstartplot+1;

figure
shadedErrorBar(tvec_new,distvec_smooth_new,distste_new,'lineProps',{'-k','lineWidth',2})
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('v(t)')
xlabel('Time (min)')
ax = gca;
ylim([4.5,6])
line([tvec_new(onset_new),tvec_new(onset_new)],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
% xlim([0,max_ck_real+epc_postasleep+1])
xlim([-60,10])
xticks([-60:10:10])
% 
% xticks([40:200:max_ck_real+epc_postasleep+1])
% xticklabels([-70:10:10])

yl = ylim;
p_th = 0.025 /length(pvec_dist);

for iii = 1:length(pvec_dist)
    tidx = iii+basetest_epc-1-idxstartplot;
    
    if pvec_dist(iii) < p_th   % Significant
        line([tvec_new(tidx),tvec_new(tidx+1)],[0.98*yl(2),0.98*yl(2)],'color','k','linewidth',1)
    end
end


%% Figure 2C

clear all
clc
load Figure2C.mat

distall_avg_smooth = smoothdata(dotp_avg,"movmean",20);
p_th = 0.025 /(length(pvec_dist)-20);

figure
hold on
shadedErrorBar(tvec,distall_avg_smooth(time_start_sampenough:end),dotp_ste(time_start_sampenough:end),'lineProps',{'-k','lineWidth',2})
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('Alignment')
xlabel('Time (min)')
ax = gca;
% ylim([0,11])
ylim([-0.04,0.1])
yticks([-0.04:0.02:0.1])
line([tvec(max_ck_real),tvec(max_ck_real)],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
% xlim([0,max_ck_real+epc_postasleep+1])
yl = ylim;

xlim([-60,10])
xticks([-60:10:10])
% xticks([62:200:max_ck_real+epc_postasleep+1])
% xticklabels([-70:10:10])

for iii = 1:length(pvec_dist)
    tidx = iii+basetest_epc-1;
    
    if pvec_dist(iii) < p_th   % Significant
        line([tvec(tidx),tvec(tidx+1)],[0.95*yl(2),0.95*yl(2)],'color','k','linewidth',2)
    end
end

%% Figure 2D
% This figure uses data from previous

clear all
clc
load Figure2C.mat

tvec = -(max_ck_real-1)/20:0.05:10.5;
dotp_all = smoothdata(dotp_avg,"movmean",10);

load Figure2D.mat
% distall_avg1 = distall_avg - min(distall_avg);
distsmooth = smoothdata(distall_avg,2,'movmean',5);
distsmooth = distsmooth -min(distsmooth);

figure
scatter(distsmooth(1:end-20),(dotp_all.^2) /2,20,'k','filled')
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
% ylim([0,7e-4])
% xlim([0,5])
% yticks(0:1e-4:7e-4)
xlabel('s')
ylabel('E (Potential)')

%% Figure 2E

% To plot the theoretical bifurcation diagram, the solutions of the
% differential equations are required; 
%
% Warning: this part of the code runs for long

clear all
clc
load Figure2A.mat

syms x
param_min = params_optim;
c = dd(:,2);
r = param_min(1);
K = param_min(2);
h= param_min(3);
m = param_min(4);

idx1 = find(c == c_3sol(1));
idx2 = find(c == c_3sol(end));

figure
hold on
% Plot in different ranges
cvec = [];
xvec = [];
for idxc = 1:idx1-1
    
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
    if length(xsol) ~= 1
        warning('Check code')
    end
    cvec = [cvec,c(idxc)];
    xvec = [xvec,xsol];

end
plot(cvec, xvec,'k','lineWidth',2);

cvec = [];
xvec = [];
for idxc = idx1:idx2
    
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
    if length(xsol) ~= 3
        warning('Check code')
    end
    cvec = [cvec,c(idxc)];
    xvec = [xvec;xsol];

end
plot(cvec, xvec(:,1),'k','lineWidth',2);
plot(cvec, xvec(:,2),'k--','lineWidth',2);
plot(cvec, xvec(:,3),'k','lineWidth',2);

cvec = [];
xvec = [];
for idxc = idx2+1:length(c)
    
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
    if length(xsol) ~= 1
        warning('Check code')
    end
    cvec = [cvec,c(idxc)];
    xvec = [xvec,xsol];

end
plot(cvec, xvec,'k','lineWidth',2);

xlabel('Control parameter')
ylabel('System state')
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
line([c_critc,c_critc],ylim,'LineStyle','--','LineWidth',2,'Color','r')

%% Figure 2F

clear all
clc
load Figure2F.mat

p_th = 0.025 / length(pvec_slp);         % Bonferroni correction
% p_th = 0.025;
tews_plt = (tews)/60;

xl = [-25,10];
xticksthis = [-25:5:10];

idxplt = find(tews_plt<10);
tews_plt = tews_plt(idxplt);

onsetmark = find(tews_plt == 0);      % Sleep onset marker

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
title('Average sleep stages')
ylabel('Sleep stages')
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylim([-0.5,3])
yl = ylim;
xlim(xl)
xticks(xticksthis)
line([tews_plt(onsetmark),tews_plt(onsetmark)],yl,'LineStyle','--','LineWidth',2,'Color','r')
xl = xlim;
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
ylim([0.4,0.65])
yticks([0.4:0.05:0.65])
yl = ylim;
ylabel('Autocorrelation')
line([tews_plt(onsetmark),tews_plt(onsetmark)],yl,'LineStyle','--','LineWidth',2,'Color','r')
pvecnow = pvec_ews{1};
for iii = 1:length(pvecnow)
    tidx = iii+base_epcs-1-lengap;
    
    if pvecnow(iii) < p_th   % Significant
        line([tews_plt(tidx),tews_plt(tidx+1)],[0.98*yl(2),0.98*yl(2)],'color','k','linewidth',1)
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
ylabel('Standard deviation')

xlabel('Time (min)')
pvecnow = pvec_ews{2};
for iii = 1:length(pvecnow)
    tidx = iii+base_epcs-1-lengap;
    
    if pvecnow(iii) < p_th   % Significant
        line([tews_plt(tidx),tews_plt(tidx+1)],[0.99*yl(2),0.99*yl(2)],'color','k','linewidth',1)
    end
end
yl = ylim;
yl = [1.2,1.8];
line([tews_plt(onsetmark),tews_plt(onsetmark)],yl,'LineStyle','--','LineWidth',2,'Color','r')

%% Figure 2G

clear all
clc
load Figure2G.mat

figure
scatter(depth_pos,r2all,20,'filled')
p = polyfit(depth_pos,r2all,1);
hold on
xl = xlim;
dx = diff(xl)/100;
xvec = xl(1):dx:xl(2);
plot(xvec,xvec*p(1)+p(2),'r','LineWidth',2)
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Avg. Sleep score')
ylabel('Fitting accuracy (R2)')
ylim([-1.5,1])
yticks([-1.5:0.5:1])

%% Figure 2H

clear all
clc
load Figure2H.mat

f = figure;
f.Position(3:4) = [500,600];

tonset = find(tvec==0);

ax1=subplot(3,1,1);
plot(tvec,hyp_expand,'LineWidth',2)
ylim([0,3])
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
set(gca, 'YTick', 0:3, 'YTickLabel', stageNames);
line([0,0],ax1.YLim,'LineStyle','--','LineWidth',2,'Color','r');

ax2=subplot(3,1,2:3);
plot(tvec,s_dynam,'k')
hold on
plot(tvec,sfit,'b','LineWidth',3)
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('s')
xlabel('Time (min)')
line([0,0],ax2.YLim,'LineStyle','--','LineWidth',2,'Color','r');


%% Figure 2I

clear all
clc
load Figure2I.mat

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
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Fitting accuracy (R2)')
ylabel('No. of participants')
xlim([0,1])
r2med = median(r2all_new);
line([r2med,r2med],ylim,'LineStyle','--','LineWidth',2,'Color','r');















