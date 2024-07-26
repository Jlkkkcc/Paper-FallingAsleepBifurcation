function [sig_artnan, ecg_artnan, sig_rm, ecg_rm,score_rm,score_rm_clean] = genart_rm_artnan(sig,ecg,Fs,epc_len,num_rep,score)

% This function serves as a general artefact rejection before the pipeline
% starts; Methods are simply using the RMS thresholding;
%
% Author: Junheng Li
%
% Inputs:
% sig: the EEG signal that needs artefact rejection; can be multi-channel
% signal NxM where N is the number of channels;
% ecg: the input ECG signal to be cleaned together
% Fs: the sampling rate of the signal
% epc_len: the epoch length basis for artefact detection in 's'
% num_rep: number of repetitions in detection
% score: the relevant scoring of the signal based on the epochs given
%
% Outputs:
% sig_rm: signals after artefact rejection with artefact epochs set to NaNs
% score_rm: the artefact-existing epoch will be set to NaN;

epc_pt = epc_len*Fs;

len_EEG = size(sig,2); 
num_epc = floor(len_EEG/epc_pt);

t = 0:1/Fs:(len_EEG-1)/Fs;

num_ch = size(sig,1);
score_rm = score(1:num_epc);
sig_rm = sig;
ecg_rm = ecg;

art_ch = zeros(num_ch,num_epc);

for ch = 1:num_ch
    
    EEG_ch = sig(ch,:);   % EEG trace of current channel;
    
    rms_all_ch = zeros(num_epc,1);
    start_epc = zeros(num_epc,1);
    
    for i = 1:num_epc
        start = 1+(i-1)*epc_pt;
        eeg_epc = EEG_ch(start:start+epc_pt-1);
        
        rms_all_ch(i) = rms(eeg_epc);
        start_epc(i) = start;
        
    end
    for ite = 1: num_rep
    
        med_rms = nanmedian(rms_all_ch);
        th = 2*med_rms;
        for i = 1:num_epc
            if rms_all_ch(i)> th
                start = 1+(i-1)*epc_pt;
                art_ch(ch,i) = 1;
                score_rm(i) = NaN;
                rms_all_ch(i) = NaN;
                sig_rm(:,start:start+epc_pt-1) = NaN(num_ch,epc_pt);    % Directly delete the artefact epochs
                ecg_rm(start:start+epc_pt-1) = NaN(1,epc_pt);
            end
        end
        
    end
   
end

sig_artnan = sig_rm;
ecg_artnan = ecg_rm;

idx_nan = isnan(sig_rm(1,:));
sig_rm(:,idx_nan) = [];
ecg_rm(isnan(ecg_rm)) = [];

score_rm_clean = score_rm;
score_rm_clean(isnan(score_rm_clean)) = [];    % Cleaned scoring





