
%% Example code for bifurcation function fitting

% Author JL

% This code gives an exemplary data for one to explore different parameters
% that can affect the function fit towards the Sleep-Distance dynamics

% You can explore around the initial conditions and see how it can affect
% fittings; Then you can use the optimisation function we provided to see
% our optimised fit

%% Data 

clear all
clc
load SleepDistance_Example.mat

% Variable explanation
% tvec_plot: the Time vector of the Sleep distance dynamics
% xx_smooth: the Sleep distance dynamics

%% Trial of parameters

% Initial parameters set-up
% You can change the parameters to see differences
K = mean(xx_smooth(1:100));
r = 3;
h = 0.6;
m = 0.45;
params = [r,K,h,m];

% Initial Condition
x_ini = [mean(xx_smooth(1:100));1];

% Dynamics from this parameter
[t,dd] = ode45(@(t,x) harvest(t,x,params),tvec_plot,x_ini);
figure
% plot(t,dd(:,1),'-o',t,dd(:,2),'-.')
plot(t,dd(:,1),'-o')
hold on
plot(t,xx_smooth,'Color','k')
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
ylabel('Sleep Distance')
xlabel('Time (min)')

% Fit R-square under this parameter ; bad fit?

rsq = rsquare(xx_smooth,dd(:,1))


%% Optimisation from the initial parameter
% Note that given our optimisation is not strictly grid search, if the
% initial parameters are extremely bad or far away from a potential fit,
% the function will fail.

% You can set the warnings to off, because the ODE solvers will excede
% tolerances very often.
warning off

% First-step tuning, tune the K values
% If you do not input the initial parameters, the function will have a
% default value
[param_tuned,rsq_init,rsq_final,dd,xini,iffail] = tunebif_param(params,xx_smooth,tvec_plot);
figure
plot(tvec_plot,dd(:,1),'-o')
hold on
plot(tvec_plot,xx_smooth)
rsq_final

% Second step, grid searching parameters
[params_optim,rsq_ini,mse_init,rsq_final,mse_final,dd,iffail1] = finetune_bif_rev(param_tuned,xx_smooth,tvec_plot,xini,[]);

figure
plot(tvec_plot,dd(:,1),'-o')
hold on
plot(tvec_plot,xx_smooth)
rsq_ini 
rsq_final

%% Tipping point evaluation

[SVariableCriticalZone_c, c_3sol_output_c] = plot_bifurcation_rev(dd, xx_smooth, tvec_plot, 'Noname', params_optim, 1);


