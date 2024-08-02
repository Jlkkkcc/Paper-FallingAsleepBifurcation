
%% This code shows a complete prediction outcome of one subject:

% Author: Anastasia Ilina / Junheng Li

% One training night, and all the remaining nigts predicted;
%%%%%%%%
% Feel free to play around with this data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable explanation:
% dd: The Bifurcation function fit output; First column is the fit of
% System state (Sleep distance dynamics), and the second column is the
% control parameter of the fit;

% estimatedSVariables: The predicted Sleep Distance dynamics of all testing
% nights;

% params_optim: Optimal bifurcaiton fit parameters of the training night

% Tcross_raw: Post-hoc model tipping point times for all testing nights

% tvec: time vector

% xx_smooth: The smoothed Sleep distance dynamics of the training night

% criticalSValue: The critical value threshold (for predicting tipping
% points) derived from the training night model fitting


%% Data

clear all
clc

load PredictionExample.mat
load Figure4D1.mat


%% Plotting training nights and prediction nights as a whole

% Controls the y-axies limit
yl = [0,12];

% Plot all testing nights
% Blue line shows the Model tipping point time; Red line shows hypnogram
% sleep onset; The crossing detections are shown in red markers;

% The methodoloy to detect crossings are embedded in this function as well
plot_predictions_together(estimatedSVariables, criticalSValue, globalMin,yl,Tcross_raw)

% Training night and its bifurcation function fit
[SVariableCriticalZone, c_3sol_output]  = plot_bifurcation_rev(dd, xx_smooth, tvec, 'Noname', params_optim, 1);


