%% Early warning signal (EWS) evaluation on continuous sleep EEG
% Version Number 3: @22/06/2022
% Modified from previous code, but now the input time series can have NaN
% valuee; The aim is to check whether median filtering is gonna destroy the
% original time-series non-linearity and cause spurious correlations.
%
% Author: Junheng Li;
% Most methods are referred from paper:
% Dakos, Vasilis, et al. "Methods for detecting early warnings of critical 
% transitions in time series illustrated using simulated ecological data." 
% PloS one 7.7 (2012): e41010.
% 
% Partially adapted from their toolbox 
% https://git.wur.nl/sparcs/generic_ews-for-matlab/-/tree/master
%
% Adapted and incorporate the analysis capability for continuous and long
% sleep time series data (original EEG, features)

%% Main codes

function [ews_out, tau] = ews_sleep_paper(ts,time,win,ifdetrend)

% The main EWS function
% Inputs:
% ts: The continuous sleep time-series data
% time: The real time index for calculating Kendall's tau, should be of the
% same length as ts
% win: The window size for EWS evaluation (Number of data points)
% ifdetrend: If to detrend the data before EWS evaluation (Linear
% detrending used here)
%
% EWS algorithms and outputs:
% ews_out: The data structure containing EWS measurements:
%          AR1 and AR2: autocorrelation at lag1 and lag2
%          Var and CoV: Variance and Coefficient of variation
%          DFA: Detrended fluctuation analysis (best polynomial fit
%          exponent)
%          TVAR: Time varying autoregression analysis; The parameter used
%          is the characteristic root of regression model;
%          Ives, A. R., and V. Dakos. 2012. Detecting dynamical changes in
%          nonlinear time series using locally linear state-space models.
%          Ecosphere 3(6):58. http://dx.doi.org/10.1890/ES11-00347.1 
% tau: The Kendall's tau of the EWS measurements (as a correlation with
% time); NO tau measurment for TVAR;


len = length(ts);
if win>=len
    error('Too short input timeseries')
end

% Lag-1 and Lag-2 autocorrelation measurements
acf1_ts = zeros(1,len-win+1);
acf1_ts_new = zeros(1,len-win+1);
% Coefficient of variation 
std_ts = zeros(1,len-win+1);
% DFA
dfa_ts = zeros(1,len-win+1);
% Time varying AR analysis
% tvar_ts = zeros(1,len-win+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EWS evaluation using rolling windows
for i = win:len

    ts_win = ts((i-win+1):i);            % Data in current window

    if ifdetrend
        ts_win = detrend(ts_win,1);    % Linear detrending
    end
    

    % Check if this window are all NaNs, if so return error
    if sum(1-isnan(ts_win)) ==0
        error('The number of continuous NaNs should be smaller than the running window size')
    end
    
    % Autocorrelation
    acf1_ts(i-win+1) = acf(ts_win,1);

    % Std
    std_ts(i-win+1) = nanstd(ts_win);
    
    % DFA
    dfa_ts(i-win+1) = dfa(ts_win);
    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EWS evaluation (model-based)
% Time varying autoregression (Characteristic root)
% AR model with order 1 used only 
% ts_log = log(ts+1);
% ts_zslog = (ts_log-nanmean(ts_log))/std(ts_log);
% options.figplot = 0;
% [lnlike, ~, ~, ~, ~, ~, ~, ~, predictvals] = TVARSS(ts_zslog, 1,options);
% tvar_ts = predictvals(:,end);          % The inverse of the roots

% All measurements
ews_out.AR1 = acf1_ts;
ews_out.StD = std_ts;
ews_out.DFA = dfa_ts;

% ews_out.TVAR = tvar_ts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kendall's tau evaluation for most metrics
% Notice that tau should be calculated against real time indexes
t_ts = time(win:end);
if size(t_ts,1)>1
    t_ts = transpose(t_ts);
end
% if size(time)>1
%     time = transpose(time);
% end

tau.AR1 = corr(t_ts', acf1_ts','Type','Kendall');
tau.StD = corr(t_ts', std_ts','Type','Kendall');
tau.DFA = corr(t_ts', dfa_ts','Type','Kendall');
% tau.TVAR = corr(time', tvar_ts,'Type','Kendall');


end

%% Autocorrelation calculation function
function res = acf(X, alag)
    if nargin < 2
        alag = 1;
    end
    %autocorrelation function
    n = length(X);
    s = var(X,'omitnan');
    mu = mean(X,'omitnan');
    Xt = X(1:end - alag);
    Xtk = X(1 + alag:end);
    res = 1 / (n - 1) / s * sum((Xt - mu) .* (Xtk - mu),'omitnan');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detrended fluctuation analysis
function res = dfa(y)         % Finding the optimal polynomial fitting exponent
    
% DFA implemented and adapted from :Livina, V. N., and Lenton, T. M.
% (2007), A modified method for detecting incipient bifurcations in a
% dynamical system, Geophys. Res. Lett., 34, L03712,
% doi:10.1029/2006GL028672.

if size(y,1)==1             % Make sure y is column vector
    y = transpose(y);
end

Yprof = cumsum(y,'omitnan');          % time series profile
lenY = length(Yprof);
smin = 10;          % minimum DFA segment length
smax = floor(lenY/2);
if lenY <= smin
    warning('Window length too short for DFA analysis, returning NaN')
    res = NaN;
    return
end
if smin>=smax
    warning('Window length too short for DFA analysis, returning NaN')
    res = NaN;
    return
end
sstep = 2;

ns = 0;
svec = [];
for s = smin:sstep:smax
    
    ns = ns +1;
    s = floor(s);
    t = (1:s)';
    svec(ns) = s;
    Nseg = floor(lenY/s);
    ordDFA = 1;
    Yerr = zeros(s,Nseg);
    for i = 1:Nseg
        Yseg = Yprof((i-1)*s+1:i*s);
        p = polyfit(t,Yseg,ordDFA);
        Yfit = polyval(p,t);
        Yerr(:,i) = Yseg-Yfit;
    end
    
    F(ns) = sqrt(mean(sum(Yerr.^2)/s));

end

pF = polyfit(log10(svec),log10(F),1);
res = pF(1);

end

