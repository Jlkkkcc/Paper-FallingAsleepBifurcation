%% Supplementary figures 6-7

% Author: JL

%% Supplementary Figure 6

clear all
clc
load SuppFig6.mat

p_th = 0.025/length(pvec_dist);
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
ylim([0,2])
xlim([-60,10])
xticks([-60:10:10])
line([tvec_now(onsetmark),tvec_now(onsetmark)],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
% xlim([0,max_ck_real+epc_postasleep+1])
yl = ylim;

scatter(tvec_now(idx_critic),dd(idx_critic,1),80,'red','filled')
% 
% xticks([62:200:max_ck_real+epc_postasleep+1])
% xticklabels([-70:10:10])

for iii = 1:length(pvec_dist)
    tidx = iii+basetest_epc-1-idxstartfit;
    
    if pvec_dist(iii) < p_th   % Significant
        line([tvec_now(tidx),tvec_now(tidx+1)],[0.95*yl(2),0.95*yl(2)],'color','k','linewidth',1)
    end
end

%% Supplementary Figure 7  

% To get the same figures as in the paper, please choose the feature
% combination to plot; We did not list them exhaustively;
%
% Revise the x and y axis limits for different features used to generate
% the right plots

% The four key feature indexes used in the paper:
% FPC1: Total EEG power (41), Spectrum centroid (46)
% FPC2: Delta-to-Alpha ratio (6), Theta temporal coherence (17)

clear all
clc
load SuppFig7.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please choose the feature combination to plot 
ftn = [6,17,41,46];     % All key FPCA feature index

% Which 2 feature combination to plot
fttest1 = 3;
fttest2 = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time to plot

epc_len = 6;
jmp_len = 3;
score_size = 30;

time_before = 60*60;   
epcs_before = floor((time_before-epc_len)/jmp_len)+1;

post_per = 21;
post_epcs = floor((post_per*score_size-epc_len)/jmp_len)+1;

maxtslen = epcs_before + post_epcs +1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ftts_avg_norm = ftall_mat_allnorm_strm;
% tvec = (-time_ews+6:dt:8*60);
tvec = t_now_onset;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Find median of each time period and evaluate direction vectors

gap_eval = 30;           % Gap of measurement
ftidx = [fttest1,fttest2];
gap_eval_dp = floor((gap_eval)/jmp_len);

all_points = [ftts_avg_norm(ftn(fttest1),end - maxtslen+1:end);ftts_avg_norm(ftn(fttest2),end - maxtslen+1:end)]';

num_steps = floor(maxtslen/gap_eval_dp);

grid_coord = [];
all_vector = [];
vec_median = [];
tvec_coord = [];

for nstp = 1:num_steps-1
    
    idxstep = (nstp-1)*gap_eval_dp+1;
    idx_nextstp = (nstp)*gap_eval_dp+1;

    cur_points = all_points(idxstep:idx_nextstp-1,:);
    nextpts = all_points(idx_nextstp:idx_nextstp+gap_eval_dp-1,:);

    med_cur = median(cur_points);
    grid_coord = [grid_coord;med_cur];
    
    % Time vector of medians
    tvec_coord = [tvec_coord;t_now_onset(idx_nextstp)-6];

    med_next = median(nextpts);

    diff_cur = nextpts - cur_points;
    for iii = 1:2
        diff_cur(:,iii) = smoothdata(diff_cur(:,iii),1,'Gaussian',10);
    end

    all_vector = [all_vector;mean(diff_cur)];

    vec_median = [vec_median;med_next - med_cur];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot median points only

figure
scatter(grid_coord(:,1),grid_coord(:,2),80,tvec_coord/60,'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
hold on
quiver(grid_coord(:,1),grid_coord(:,2),vec_median(:,1),vec_median(:,2),...
    'color',[0,0,0],'linewidth',1.5)

idxonsethere = find(tvec_coord==0);
pltonsetmark = idxonsethere:idxonsethere+1;

% Add special color to sleep onset
scatter(grid_coord(pltonsetmark,1),grid_coord(pltonsetmark,2),80,'red','filled')
gg = parula(length(tvec_coord));
gg(pltonsetmark,:) = repmat([1,0,0],length(pltonsetmark),1);
colormap(gg)
c = colorbar;
c.Limits = [-60,10];
c.Ticks = [-60:10:10];
c.Label.String = 'Time (min)';

set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel(ft_dscrp{ftn(fttest1)})
ylabel(ft_dscrp{ftn(fttest2)})


% Control y-axis lim
% For different feature combinations this need to be revised
%%%%%%%%%%%%%%%%
% ylim([-1.8,0.2])
% yticks([-1.8:0.4:0.2])
% 
% % ylim([0,0.7])
% 
% ylim([-0.2,1.8])
% yticks([-0.2:0.4:1.8])
% xlim([0,0.7])
xlim([-1.8,0.2])
xticks([-1.8:0.4:0.2])










