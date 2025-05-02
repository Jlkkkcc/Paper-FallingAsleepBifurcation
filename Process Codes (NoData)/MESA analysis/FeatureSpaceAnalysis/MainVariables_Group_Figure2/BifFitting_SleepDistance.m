
%% Sleep distance dynamics - Group level - Bifurcation function fitting

% Author: JL

%%  Data loading

clear all
clc
load SleepDistance_Grp.mat

%% Fitting initialisation

time_start_sampenough = find(num_samps_time>=N_min,1);
max_ck_real = length(distall_avg) - time_start_sampenough - epc_postasleep;
tvec = -(max_ck_real-1)/20:0.05:10.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cut out post-sleep onset period for fitting
xx = distall_avg(time_start_sampenough:end);
xx = xx-min(xx);   % Gap correction

% % Only uses pre-asleep to fit
% % xx = xx(1:end-epc_postasleep);
% % tvec = tvec(1:end-epc_postasleep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting starting point
timefit = 60*60;
epcs_tofit = floor((timefit-epc_len)/jmp_len)+1;

idxstartfit = max_ck_real-epcs_tofit+1;

% idxstartfit = 1;

xx = xx(idxstartfit:end);
tvec = tvec(idxstartfit:end);

% Initial parameters (guess)
r = 8;
K = mean(xx(1:100))+0.5;
h = 0.5;
m = 0.089;
params = [r,K,h,m];

[t,dd] = ode45(@(t,x) harvest(t,x,params),tvec,[xx(1);1]);
figure
% plot(t,dd(:,1),'-o',t,dd(:,2),'-.')
plot(t,dd(:,1),'-o')
hold on
plot(t,xx)

%% Optimisation of fit

xx_smooth = smoothdata(xx,2,'movmedian',[20,0]);
xx_smooth = xx_smooth- min(xx_smooth);
tvec_now = tvec;

[param_tuned,rsq_init,rsq_final,dd,xini,iffail] = tunebif_param(params,xx_smooth,tvec_now);
figure
plot(tvec_now,dd(:,1),'-o')
hold on
plot(tvec_now,xx_smooth)
rsq_final

% Further fine tune

[params_optim,rsq_ini,mse_ini,rsq_final,mse_final,dd,iffail1] = finetune_bif_rev(param_tuned,xx_smooth,tvec_now,xini,[]);

figure
plot(tvec_now,dd(:,1),'-o')
hold on
plot(tvec_now,xx_smooth)
rsq_ini
rsq_final

%% Bifurcation diagram under this parameter

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

%% Critical tipping point identification

idx_critic = find(c == c_3sol(end));   % Critical point
t_critc = tvec_now(idx_critic)
c_critc = c(idx_critic)


%% Finding features (group-level) resembling to control parameter

load FeatureSpace_GroupAnalysis.mat

% Control parameter time-series
ts_range = 1:length(dd(:,2))-epc_postasleep;
% ts_range = 1:length(dd(:,2));
c_ts = dd(ts_range,2)';

c_ts_norm = c_ts./sqrt(c_ts*c_ts');

num_fttest = 47;

% Two different metrics: DTW and inner similarity
dtw_allfts = zeros(num_fttest,1);
dtw_allfts_normed = zeros(num_fttest,1);
innprod_allfts = zeros(num_fttest,1);

for iff = 1:num_fttest
    
    xnow = ftall_mat_allnorm_strm(iff,idxstartfit:end-epc_postasleep);
    x_norm = xnow./sqrt(xnow*xnow');
    
    dtw_allfts(iff) = dtw(xnow,c_ts);
    dtw_allfts_normed(iff) = dtw(x_norm,c_ts_norm);

    innprod_allfts(iff) = x_norm * c_ts_norm';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find most similar one

disp(['Minimum DTW distance: ',num2str(find(dtw_allfts == min(dtw_allfts)))])
disp(['Minimum DTW distance normed: ',num2str(find(dtw_allfts_normed == min(dtw_allfts_normed)))])
disp(['Maximum similarity: ',num2str(find(innprod_allfts == max(innprod_allfts)))])

innprod_allfts = abs(innprod_allfts);

newmet = abs(innprod_allfts)./dtw_allfts_normed;
disp(['Best New metric: ',num2str(find(newmet == max(newmet)))])

%% Save for plots

save('Figure2a.mat','epcs_tofit','distall_ste','tvec_now','dd','basetest_epc','idxstartfit','pvec_dist','params_optim','time_start_sampenough','xx_smooth','t_critc','c_critc',"c_3sol",'xx')
















