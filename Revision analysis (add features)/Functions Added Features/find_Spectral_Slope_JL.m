function [results] = find_Spectral_Slope_JL(y,Fs,f_range,python_directory,iflinux)

[psd_y,f] = pmtm(y,4,length(y),Fs);      %Multi-taper power-spectrum density analysis
psd_y_norm = psd_y ./ sum(psd_y);

exponent = extract_aperiodic_component(psd_y_norm, f, f_range, 0, python_directory, iflinux); % f_range should be in the format of [min_freq, max_freq]
results = exponent;



