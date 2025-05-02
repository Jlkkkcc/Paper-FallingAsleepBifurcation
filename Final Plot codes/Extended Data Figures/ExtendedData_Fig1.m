%% Extended Data Figures plots Figure 1

% Author: Junheng Li

%% Extended Data Fig.1a

clear all
clc
load ExtData_Fig1a.mat

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
ylim([-1,1])
yticks([-1.5:0.5:1])
xlim([0,3])
xticks([0:3])

%% Extended Data Fig.1b

clear all
clc
load ExtData_Fig1b.mat

[rr1,pp1] = corr(K_all_clean,r2clean)

figure
scatter(K_all_clean,r2clean,20,'filled')
p = polyfit(K_all_clean,r2clean,1);
hold on
xl = xlim;
dx = diff(xl)/100;
xvec = xl(1):dx:xl(2);
plot(xvec,xvec*p(1)+p(2),'r','LineWidth',2)
box off
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Distance')
ylabel('R-squared')
ylim([-1,1])
yticks([-1:0.5:1])


%% Extended Data Fig.1c

clear all
clc
load ExtData_Fig1c.mat

mean(t_crtic_new)
std(t_crtic_new)
median(t_crtic_new)

figure
edgesbin = [-10:1:10];
pd = fitdist(t_crtic_new,'Kernel');
pdfbin = 0:0.01:1;
ypdf = pdf(pd,pdfbin);

histogram(t_crtic_new,edgesbin)
% histfit(tcrtc_new,30,'kernel')
hold on
% plot(pdfbin,ypdf*(0.05*length(t_crtic_new)),'LineWidth',2,'Color','r')
box off
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Tipping point time (min)')
ylabel('No. of participants')
% xlim([0,1])
r2med = median(t_crtic_new);
line([r2med,r2med],ylim,'LineStyle','--','LineWidth',2,'Color',[0.25,0.187,0.02]);
line([0,0],ylim,'LineStyle','--','LineWidth',2,'Color','r');



%% Extended Data Fig.1d

clear all
clc
load ExtData_Fig1def.mat

[rr3,pp3] = corr((t_crtic_clean),lat_clean)

figure
scatter(lat_clean,t_crtic_clean,20,'filled')
p = polyfit(lat_clean,t_crtic_clean,1);
hold on
xl = xlim;
dx = diff(xl)/100;
xvec = xl(1):dx:xl(2);
plot(xvec,xvec*p(1)+p(2),'r','LineWidth',2)
box off
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Sleep latency (min)')
ylabel('Tipping point time (min)')
xlim([0,90])
xticks([0:10:90])
ylim([-20,10])
line(xl,[0,0],'LineStyle','--','LineWidth',2,'Color','y')


%% Extended Data Fig.1e

clear all
clc
load ExtData_Fig1def.mat

[rr2,pp2] = corr(K_all_clean,lat_clean)

figure
scatter(lat_clean,K_all_clean,20,'filled')
p = polyfit(lat_clean,K_all_clean,1);
hold on
xl = xlim;
dx = diff(xl)/100;
xvec = xl(1):dx:xl(2);
plot(xvec,xvec*p(1)+p(2),'r','LineWidth',2)
box off
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Sleep latency (min)')
ylabel('Distance')
xlim([0,90])
xticks([0:10:90])


%% Extended Data Fig.1f

clear all
clc
load ExtData_Fig1def.mat

[rr4,pp4] = corr(r2clean,lat_clean)

figure
scatter(lat_clean,r2clean,20,'filled')
p = polyfit(lat_clean,r2clean,1);
hold on
xl = xlim;
dx = diff(xl)/100;
xvec = xl(1):dx:xl(2);
plot(xvec,xvec*p(1)+p(2),'r','LineWidth',2)
box off
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Sleep latency (min)')
ylabel('Fitting accuracy (R2)')
ylim([0,1])
% yticks([0:0.5:1])
xlim([0,90])
xticks([0:10:90])


