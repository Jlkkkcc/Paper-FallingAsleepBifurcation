function [ft, ft_dscrp] = sleep_ft_extract_cont(y, Fs)

% This function is used to extract pre-defined feature set of input time
% series data;
%
% Author: Junheng Li
%
% Input:
% y: input time-series
% Fs: sampling rate
%
% Output:
% ft: feature vector
% ft_dscrp: relevant feature description

%% Log of code
% All previous, features extracted: Frequency-band features, temporal
% coherence, phase-amplitude coupling and Catch-22;

% 17/07/2020, adding peak frequency feature in each frequency band;

% 05/03/2021, adding MA filter for 1/f noise removal;

% 15/03/2021, moved the MA filter before the band-pass filter

% 05/08/2021, changing the algorithm in Temporal coherence to take the
% 'floor' of the maximum number of divisions (Instead of round)

% 05/04/2023, Code modified in that it can take inputs; Modified in that
% one can choose whether to apply 1/f noise removal; And 1/f noise removal
% won't be applied to Catch-22 features anyway;

% 24/04/2023, Code modified again that Temporal coherence evaluation is
% enhanced;

% 19/06/2023: Code finalised for this version to keep consistent with all
% analysis: 1/f noise cancellation applied to all (including CATCH22); Temporal coherence
% epochs to be 1s.

%% Start

if size(y,2)<size(y,1)     %make sure it's a row vector
    y = transpose(y);
end

%% PSD estimation part

delta_bd = [1,4];
theta_bd = [4,8];
alpha_bd = [8,12];
beta_bd = [12,30];

ft_dscrp = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency band features

ifproc = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1/f differentiation filter
y_MA = y;
clear y
y = y_MA(2:end) - y_MA(1:end-1);
y = [y,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Band pass filter

if ifproc     % Normal pre-processing
    % filter EEG, broad BP
    order = 2;
    lowFreq0 = 0.1;    % 0.1
    highFreq0 = 32;
    [b,a] = butter(order, [lowFreq0 highFreq0]/(Fs/2), 'bandpass');
    y_proc = filtfilt(b,a,y);
     
else
    y_proc = y;
end


[psd_y,f] = pmtm(y_proc,4,length(y_proc),Fs);      %Multi-taper power-spectrum density analysis

psd_y_norm = psd_y ./ sum(psd_y);    % Normalize to total power

delta_idx = (f>=delta_bd(1) & f<delta_bd(2));
theta_idx = (f>=theta_bd(1) & f<theta_bd(2));
alpha_idx = (f>=alpha_bd(1) & f<alpha_bd(2));
beta_idx = (f>=beta_bd(1) & f<beta_bd(2));

delta_power = nanmean(psd_y_norm(delta_idx));
theta_power = nanmean(psd_y_norm(theta_idx));
alpha_power = nanmean(psd_y_norm(alpha_idx));
beta_power = nanmean(psd_y_norm(beta_idx));

ft.delta = log(delta_power);
ft.theta = log(theta_power);
ft.alpha = log(alpha_power);
ft.beta = log(beta_power);

ft.DTratio = log(delta_power/theta_power);
ft.DAratio = log(delta_power/alpha_power);
ft.DBratio = log(delta_power/beta_power);
ft.TAratio = log(theta_power/alpha_power);
ft.TBratio = log(theta_power/beta_power);
ft.ABratio = log(alpha_power/beta_power);

ft.Transratio = log((delta_power+theta_power)/(alpha_power+beta_power));    % Transition power ratio

ft_dscrp = [ft_dscrp,{'delta','theta','alpha','beta','DTratio','DAratio','DBratio','TAratio','TBratio','ABratio','Transratio'}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peak Oscillation frequency
f_delta = f(delta_idx);
ft.peakdelta = f_delta(find(psd_y_norm(delta_idx) == max(psd_y_norm(delta_idx)),1));
f_theta = f(theta_idx);
ft.peaktheta = f_theta(find(psd_y_norm(theta_idx) == max(psd_y_norm(theta_idx)),1));
f_alpha = f(alpha_idx);
ft.peakalpha = f_alpha(find(psd_y_norm(alpha_idx) == max(psd_y_norm(alpha_idx)),1));
f_beta = f(beta_idx);
ft.peakbeta = f_beta(find(psd_y_norm(beta_idx) == max(psd_y_norm(beta_idx)),1));

ft_dscrp = [ft_dscrp,{'deltapeak','thetapeak','alphapeak','betapeak'}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal coherence features

epc_len = 1;     % epoch length for feature extraction in 's'

len_y = floor(length(y)/Fs);               % Modification 05/08/2021
num_epc_unit = len_y/epc_len;

for i = 1:num_epc_unit
    y_epc{i} = y((i-1)*Fs*epc_len+1:i*Fs*epc_len);
end
cnt = 0;
for i = 1:num_epc_unit-1
    for ii = i+1:num_epc_unit
        cnt= cnt+1;
        data1 = y_epc{i};
        data2 = y_epc{ii};
        [Cxy_y(cnt,:),f_cohere] = mscohere(data1,data2,[],[],0:0.5:30,Fs);
    end
end

Cxy_y_avg = nanmean(Cxy_y,1);

delta_idx = (f_cohere>=delta_bd(1) & f_cohere<delta_bd(2));
theta_idx = (f_cohere>=theta_bd(1) & f_cohere<theta_bd(2));
alpha_idx = (f_cohere>=alpha_bd(1) & f_cohere<alpha_bd(2));
beta_idx = (f_cohere>=beta_bd(1) & f_cohere<beta_bd(2));

% Temporal coherence in each frequency band
ft.TCdelta = nanmean(Cxy_y_avg(delta_idx));
ft.TCtheta = nanmean(Cxy_y_avg(theta_idx));
ft.TCalpha = nanmean(Cxy_y_avg(alpha_idx));
ft.TCbeta = nanmean(Cxy_y_avg(beta_idx));

ft_dscrp = [ft_dscrp,{'TCdelta','TCtheta','TCalpha','TCbeta'}];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase-amplitude coupling features

N = 100;
% Filter chosen is zero-phase, even order FIR filter;

d_delta = designfilt('bandpassfir','FilterOrder',N,'CutoffFrequency1',delta_bd(1),'CutoffFrequency2',delta_bd(2),'SampleRate',Fs);
y_filt_delta = filtfilt(d_delta,y);

d_theta = designfilt('bandpassfir','FilterOrder',N,'CutoffFrequency1',theta_bd(1),'CutoffFrequency2',theta_bd(2),'SampleRate',Fs);
y_filt_theta = filtfilt(d_theta,y);

d_alpha = designfilt('bandpassfir','FilterOrder',N,'CutoffFrequency1',alpha_bd(1),'CutoffFrequency2',alpha_bd(2),'SampleRate',Fs);
y_filt_alpha = filtfilt(d_alpha,y);

d_beta = designfilt('bandpassfir','FilterOrder',N,'CutoffFrequency1',beta_bd(1),'CutoffFrequency2',beta_bd(2),'SampleRate',Fs);
y_filt_beta = filtfilt(d_beta,y);

y_delta_analytic = hilbert(y_filt_delta);
y_theta_analytic = hilbert(y_filt_theta);
y_alpha_analytic = hilbert(y_filt_alpha);
y_beta_analytic = hilbert(y_filt_beta);

% Amplitudes
amp_theta = abs(y_theta_analytic);
amp_alpha = abs(y_alpha_analytic);
amp_beta = abs(y_beta_analytic);

% Phases
phase_delta = angle(y_delta_analytic);
phase_theta = angle(y_theta_analytic);
phase_alpha = angle(y_alpha_analytic);

idx_cycle_delta = find_cycles(phase_delta);
idx_cycle_theta = find_cycles(phase_theta);
idx_cycle_alpha = find_cycles(phase_alpha);

zt_DT = amp_theta(1:idx_cycle_delta).*exp(j*phase_delta(1:idx_cycle_delta));
zt_DA = amp_alpha(1:idx_cycle_delta).*exp(j*phase_delta(1:idx_cycle_delta));
zt_DB = amp_beta(1:idx_cycle_delta).*exp(j*phase_delta(1:idx_cycle_delta));
zt_TA = amp_alpha(1:idx_cycle_theta).*exp(j*phase_theta(1:idx_cycle_theta));
zt_TB = amp_beta(1:idx_cycle_theta).*exp(j*phase_theta(1:idx_cycle_theta));
zt_AB = amp_beta(1:idx_cycle_alpha).*exp(j*phase_alpha(1:idx_cycle_alpha));

ft.DTcouple = abs(sum(zt_DT)/sum(amp_theta(1:idx_cycle_delta)));
ft.DAcouple = abs(sum(zt_DA)/sum(amp_alpha(1:idx_cycle_delta)));
ft.DBcouple = abs(sum(zt_DB)/sum(amp_beta(1:idx_cycle_delta)));
ft.TAcouple = abs(sum(zt_TA)/sum(amp_alpha(1:idx_cycle_theta)));
ft.TBcouple = abs(sum(zt_TB)/sum(amp_beta(1:idx_cycle_theta)));
ft.ABcouple = abs(sum(zt_AB)/sum(amp_beta(1:idx_cycle_alpha)));

ft_dscrp = [ft_dscrp,{'DTcouple','DAcouple','DBcouple','TAcouple','TBcouple','ABcouple'}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Catch-22 cannonical features

if size(y,2)>size(y,1)     %make sure it's a column vector
    y = transpose(y);
end

[ft_values,ft_name] = catch22_all(y);
for i = 1:length(ft_values)
    ft.(ft_name{i}) = ft_values(i);
end
ft_dscrp = [ft_dscrp,ft_name];



function idx = find_cycles(phase)
stp = phase(1);
for i = length(phase):-1:2
    if phase(2)>stp
        if phase(i) >= stp && phase(i-1)<= stp
            idx = i;
            return
        end
    end
    if phase(2)<stp
        if phase(i) <= stp && phase(i-1)>= stp
            idx = i;
            return
        end
    end
end
