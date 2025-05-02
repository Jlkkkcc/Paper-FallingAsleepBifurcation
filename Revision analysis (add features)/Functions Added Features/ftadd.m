
%% Additional feature evaluation function

% JL

% Revisions for paper December 2024

% Adding three features: Sigma band power / Spectral Slope / Normalised LZ
% complexity

%%

function [ft,ft_dscrp] = ftadd(y,Fs)


%% Start

if size(y,2)<size(y,1)     %make sure it's a row vector
    y = transpose(y);
end

%% PSD estimation part

sigma_bd = [12,16];

ft_dscrp = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency band features

ifproc = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1/f differentiation filter
y_MA = y;
y_1f = y_MA(2:end) - y_MA(1:end-1);
y_1f = [y_1f,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Band pass filter

if ifproc     % Normal pre-processing
    % filter EEG, broad BP
    order = 2;
    lowFreq0 = 0.1;    % 0.1
    highFreq0 = 32;
    [b,a] = butter(order, [lowFreq0 highFreq0]/(Fs/2), 'bandpass');
    y_proc = filtfilt(b,a,y_1f);
     
else
    y_proc = y_1f;
end


[psd_y,f] = pmtm(y_proc,4,length(y_proc),Fs);      %Multi-taper power-spectrum density analysis

psd_y_norm = psd_y ./ sum(psd_y);    % Normalize to total power

sigma_idx = (f>=sigma_bd(1) & f<sigma_bd(2));

sigma_power = mean(psd_y_norm(sigma_idx),'omitmissing');


ft.sigma = log(sigma_power);

ft_dscrp = [ft_dscrp,{'delta_sum','theta_sum','alpha_sum','beta_sum','sigma'}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LZ Complexity

binarized_y = y > mean(y,2);
lempel_ziv = LZ76(binarized_y);
entropy_rate_lz = lempel_ziv*log2(length(binarized_y))/length(binarized_y);

ft.LZcomp = entropy_rate_lz;
ft_dscrp  = [ft_dscrp,'LZComplexity'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1/f slope

py = pyenv;
exponent = find_Spectral_Slope_JL(y,Fs,[0.1 32],py.Executable,1);

ft.Exp1f = exponent;
ft_dscrp  = [ft_dscrp,'Exp1f'];

%%%%%%%%%%%%%%%%%%%%




end









