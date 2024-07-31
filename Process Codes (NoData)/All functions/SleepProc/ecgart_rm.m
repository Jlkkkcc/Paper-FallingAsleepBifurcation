function [sig_rm, avg_EEG_sign_out, avg_EEG_sign_res_out] = ecgart_rm(sig,ecg,Fs,RR_th)

% This function is a realization of ECG artefact removal with methods
% described in paper:
% Purcell, S., Manoach, D., Demanuele, C. et al. Characterizing sleep
% spindles in 11,630 individuals from the National Sleep Research Resource.
% Nat Commun 8, 15930 (2017). https://doi.org/10.1038/ncomms15930
%
% Author: Junheng Li
%
% Input:
% sig: the input EEG signal for ECG artefact removal; Can be a NxM matrix
% where N is the number of channels
% ecg: the ECG signal trace
% Fs: the sampling frequency
% RR_th: the amplitude threshold of detecting RR peaks in ECG
% 
% Output
% sig_rm: the clean signal without ECG artefacts
% avg_EEG_sign: original EEG trace
% avg_EEG_sign_res: the residual ECG artefact trace after removal

[~, R_t, ~, R_index, ~, ~]  = rpeakdetect(ecg',Fs,RR_th,0);

num_ch = size(sig,1);
sig_rm = sig;
len_max = 2*Fs;    % maximum signature length

avg_EEG_sign_out = zeros(num_ch,len_max);

avg_EEG_sign_res_out = zeros(num_ch,len_max);

for ch = 1:num_ch
    
    eeg_now = sig(ch,:);

    % Calculating signatures
    
    forward_shift = 10;              % Number of data points forwarded in each time-locked event
    EEG_sign = zeros(length(R_t)-1,len_max);
    
    for i = 1:length(R_t)-1
        intv = R_index(i+1) - R_index(i);
        if R_index(i) <= forward_shift
            EEG_sign(i,:) =  NaN(1,len_max);
            continue;
        end
        if intv >= len_max
            EEG_sign(i,:) = eeg_now(R_index(i)-forward_shift:(R_index(i) + len_max-forward_shift-1));
        else
            EEG_sign(i,1:intv) = eeg_now(R_index(i)-forward_shift:(R_index(i) + intv-forward_shift-1));
        end
        
    end
    
    avg_EEG_sign = nanmean(EEG_sign);
    avg_EEG_sign_out(ch,:) = avg_EEG_sign;
    
    % Average sign removal
    
    for i = 1:length(R_t)-2
        R_int = R_index(i+1) - R_index(i);
        if R_int>len_max
            R_int = len_max;
        end
        if R_index(i) <= forward_shift
            continue;
        end
        sig_rm(ch,(R_index(i)-forward_shift):(R_index(i) + R_int-forward_shift-1)) = eeg_now((R_index(i)-forward_shift):(R_index(i) + R_int-forward_shift-1)) - avg_EEG_sign(1:R_int);
    end
    
    % Residual artefact quantification
    
    EEG_sign_res = zeros(length(R_t)-1,len_max);
    eeg_new = sig_rm(ch,:);
    
    for i = 1:length(R_t)-1
        intv = R_index(i+1) - R_index(i);
        if R_index(i) <= forward_shift
            EEG_sign_res(i,:) = NaN(1,len_max);
            continue;
        end
        if intv >= len_max
            EEG_sign_res(i,:) = eeg_new(R_index(i)-forward_shift:(R_index(i) + len_max-forward_shift-1));
        else
            EEG_sign_res(i,1:intv) = eeg_new(R_index(i)-forward_shift:(R_index(i) + intv-forward_shift-1));
        end
        
    end
    
    avg_EEG_sign_res = nanmean(EEG_sign_res);
    avg_EEG_sign_res_out(ch,:) = avg_EEG_sign_res;
    
end

