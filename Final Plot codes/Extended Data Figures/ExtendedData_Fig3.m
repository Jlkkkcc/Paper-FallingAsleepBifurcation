%% Extended Data Fig3

% JL

%% Extended Data Fig3a

clear all
clc
load ExtData_Fig3a.mat

p_th = 0.025 /length(pvec_dist);

tvec = -(max_ck_real-1)/20:0.05:10.5;

figure
shadedErrorBar(tvec,distall_avg_smooth(time_start_sampenough:end),distvec_ste_newalign(time_start_sampenough:end),'lineProps',{'-k','lineWidth',2})
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('Distance')
xlabel('Time (min)')
ax = gca;
ylim([0,5])
line([tvec(max_ck_real),tvec(max_ck_real)],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
% xlim([0,max_ck_real+epc_postasleep+1])

yl = ylim;
xlim([-30,10])

% xticks([1:100:max_ck_real+epc_postasleep+1])
% xticklabels([-35:5:10])

for iii = 1:length(pvec_dist)
    tidx = iii+basetest_epc-1;
    
    if pvec_dist(iii) < p_th   % Significant
        line([tvec(tidx),tvec(tidx+1)],[0.95*yl(2),0.95*yl(2)],'color','k','linewidth',1)
    end
end

%% Extended Data Fig3b

clear all
clc
load ExtData_Fig3b.mat

p_th = 0.025 /length(pvec_dist);

tvec = -(max_ck_real-1)/20:0.05:10.5;

% distall_avg_smooth = smoothdata(distvec_avg_newalign,2,'movmedian',5);
figure
shadedErrorBar(tvec,distall_avg_smooth(time_start_sampenough:end),distvec_ste_newalign(time_start_sampenough:end),'lineProps',{'-k','lineWidth',2})
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('State velocity')
xlabel('Time (min)')
ax = gca;
ylim([5,6.5])
line([tvec(max_ck_real),tvec(max_ck_real)],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
% xlim([0,max_ck_real+epc_postasleep+1])

yl = ylim;
xlim([-30,10])

% xticks([1:100:max_ck_real+epc_postasleep+1])
% xticklabels([-35:5:10])

for iii = 1:length(pvec_dist)
    tidx = iii+basetest_epc-1;
    
    if pvec_dist(iii) < p_th   % Significant
        line([tvec(tidx),tvec(tidx+1)],[0.98*yl(2),0.98*yl(2)],'color','k','linewidth',1)
    end
end


%% Extended Data Fig3c

clear all
clc
load ExtData_Fig3c.mat

tvec = -(max_ck_real-1)/20:0.05:10.5;
p_th = 0.025 /length(pvec_dist);
figure
shadedErrorBar(tvec,distall_avg_smooth(time_start_sampenough:end),distvec_ste_newalign(time_start_sampenough:end),'lineProps',{'-k','lineWidth',2})
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('Sleep Velocity')
xlabel('Time (min)')
ax = gca;
% ylim([4,7])
% ylim([0,8])
ylim([-0.05,0.1])
line([tvec(max_ck_real),tvec(max_ck_real)],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
xlim([-30,10])

% xticks([98:100:max_ck_real+epc_postasleep+1])
% xticklabels([-30:5:10])

yl = ylim;


for iii = 1:length(pvec_dist)
    tidx = iii+basetest_epc-1;
    
    if pvec_dist(iii) < p_th   % Significant
        line([tvec(tidx),tvec(tidx+1)],[0.95*yl(2),0.95*yl(2)],'color','k','linewidth',1)
    end
end

%% Extended Data Fig3d

clear all
clc
load ExtData_Fig3d.mat

figure
shadedErrorBar(1:max_ck_real+epc_postasleep+1,distall_avg(time_start_sampenough:end)-minall,distall_ste(time_start_sampenough:end),'lineProps',{'-k','lineWidth',2,'Color',[0,0,0,1]})
hold on
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('Sleep Distance')
xlabel('Time (min)')
ax = gca;
ylim([0,5])
line([max_ck_real,max_ck_real],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
xlim([max_ck_real-1200,max_ck_real+epc_postasleep+1])

iddst = max_ck_real-1200;
shadedErrorBar(max_ck_real-max_ck_real11+1:max_ck_real+epc_postasleep+1,distall_avg11(time_start_sampenough11:end)-minall,distall_ste11(time_start_sampenough11:end),'lineProps',{'-b','lineWidth',2,'Color',[0,0,1,0.2]})

xticks([iddst:200:max_ck_real+epc_postasleep+1])
xticklabels([-60:10:10.5])


yl = ylim;
% p_th = 0.05;
p_th = 0.05 / length(time_start_sampenough11:max_epctotal);
for iii = time_start_sampenough11:max_epctotal
    tidx = max_ck_real-max_ck_real11+1 + iii - time_start_sampenough11;
    
    if pval(iii) < p_th   % Significant
        line([tidx,tidx+1],[0.98*yl(2),0.98*yl(2)],'color','k','linewidth',1)
    end
end
% 

%% Extended Data Fig3e

clear all
clc
load ExtData_Fig3e.mat

figure
shadedErrorBar(1:max_ck_real+epc_postasleep+1,distall_avg_smooth(time_start_sampenough:end),distall_ste(time_start_sampenough:end),'lineProps',{'-k','lineWidth',2})
hold on
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('State velocity')
xlabel('Time (min)')
ax = gca;
ylim([4.5,6])
line([max_ck_real,max_ck_real],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
xlim([max_ck_real-1200,max_ck_real+epc_postasleep+1])

iddst = max_ck_real-1200;
shadedErrorBar(max_ck_real-max_ck_real11+1:max_ck_real+epc_postasleep+1,distall_avg_smooth11(time_start_sampenough11:end),distall_ste11(time_start_sampenough11:end),'lineProps',{'-b','lineWidth',2,'Color',[0,0,1,0.3]})

yl = ylim;
p_th = 0.05 / length(time_start_sampenough11:max_epctotal);
for iii = time_start_sampenough11:max_epctotal
    tidx = max_ck_real-max_ck_real11+1 + iii - time_start_sampenough11;
    
    if pval(iii) < p_th   % Significant
        line([tidx,tidx+1],[0.98*yl(2),0.98*yl(2)],'color','k','linewidth',1)
    end
end

xticks([iddst:200:max_ck_real+epc_postasleep+1])
xticklabels([-60:10:10.5])

%% Extended Data Fig3f

clear all
clc
load ExtData_Fig3f.mat

f = figure;
f.Position(3:4) = [500,500];

tonset = find(tvec==0);

ax1=subplot(3,1,1);
plot(tvec,hyp_expand,'LineWidth',2,'Color','k')
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

% Tipping point
scatter(t_crtc,sfit(find(tvec==t_crtc),1),80,'filled','MarkerFaceColor','r')






















