function exponent = extract_aperiodic_component(psd, freqs, freq_range, ifpythonreset, python_directory, iflinux)

% This function is used to extract the parameters of aperiodic component
% fit to the EEG signal usign fooof python library 
%
%
% Author: Anastasia Ilina
%
% Input:
% y: input time-series (EEG signal)
% Fs: sampling rate
% psd: normalised power spectrum for the time-series in hand (row vector)
% freqs: row vector of frequency values in the power spectrum
% f_range: fitting range (Hz)

%
% Output:
% exponent: the slope of the aperiodic component fit of the EEG signal 

%% Log of code

% 22/05/2023 - created the code 


%% Setting up Python 
%python_directory = 'usr/bin/python3';
if iflinux
    flag = int32(bitor(2, 8));
    py.sys.setdlopenflags(flag);
end 

if ifpythonreset
    pe = pyenv('Version', python_directory, "ExecutionMode","InProcess");
end 

py.importlib.import_module('numpy');
py.importlib.import_module('scipy');
py.importlib.import_module('fooof');


%% Getting the aperiodict component  
settings = struct;
fooof_results = fooof(freqs, psd, freq_range, settings, 1);
% fooof_plot(fooof_results)
exponent = fooof_results.aperiodic_params(2);