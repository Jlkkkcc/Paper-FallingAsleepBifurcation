%% Extended Data Fig7

% JL

%% Extended Data Fig7a

clear all
clc
load ExtData_Fig7a.mat

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


%% Extended Data Fig7a

clear all
clc
load ExtData_Fig7a.mat

% Warning: evaluation for bifurcation diagram takes very long

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

%% Extended Data Fig7c

clear all
clc
load ExtData_Fig7c.mat

figure
edgesbin = [0:0.05:1];
pd = fitdist(r2all(mask_all),'Kernel');
pdfbin = 0:0.01:1;
ypdf = pdf(pd,pdfbin);

histogram(r2all(mask_all),edgesbin)
% histfit(tcrtc_new,30,'kernel')
hold on
% plot(pdfbin,ypdf*(0.05*length(r2all(mask_all))),'LineWidth',2,'Color','r')
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


%% Extended Data Fig7c

clear all
clc
load ExtData_Fig7d.mat

figure
edgesbin = [-20:1:0];
pd = fitdist(TC_Mdl(mask_all),'Kernel');
pdfbin = 0:0.01:1;
ypdf = pdf(pd,pdfbin);

histogram(TC_Mdl(mask_all),edgesbin)
% histfit(tcrtc_new,30,'kernel')
hold on
% plot(pdfbin,ypdf*(0.05*length(TC_Mdl(mask_all))),'LineWidth',2,'Color','r')
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('Tipping point time (min)')
ylabel('No. of nights')
xlim([-20,0])
r2med = median(TC_Mdl(mask_all),'omitmissing');
line([r2med,r2med],ylim,'LineStyle','--','LineWidth',2,'Color','r');




