function [art_all, art_ch, proportion_epochs_remaining] = RMS_artefact_rejection(EEG_table_timed, num_epc, epc_len, stepSec, Fs, iter, th_coeff, method)

% This function is used for general artefact rejection based on large
% epochs for removing extraordinary artefacts such as movements;
%
% Author: Anastasia Ilina, Junheng Li
%
% method: either 'RMS' or 'Otsu' for selecting the thresholding technique


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-calculation
num_eeg_ch = length(EEG_table_timed.Properties.VariableNames);
eegChannels = EEG_table_timed.Properties.VariableNames; % Excluding time column
art_all = zeros(1,num_epc);
art_ch = zeros(num_eeg_ch,num_epc);

% Artefact detection
num_epcs_rejected = 0;

rms_all = zeros(num_eeg_ch,num_epc);

for j = 1:num_epc
    % Calculate the starting index for this epoch in terms of rows in the 'data' table
    startTimeIdx = (j - 1) * stepSec + 1;
    % Calculate the ending index for this epoch in terms of rows in the 'data' table
    endTimeIdx = startTimeIdx + epc_len - 1;
    
    % Debug: Print the current epoch number
    %disp(['Processing epoch: ', num2str(j)]);

    % Make sure the indices don't exceed the data limits
    if endTimeIdx <= height(EEG_table_timed)
        % Initialize a matrix to hold the EEG data for this epoch
        epochDataFallingAsleep = zeros(epc_len * Fs, length(eegChannels));
        % Extract the data for each channel for this epoch
        for ch = 1:length(eegChannels)
            channelDataFallingAsleep = [];
            

            % Extract the EEG data from consecutive cells in the table and concatenate them
            for sec = startTimeIdx:endTimeIdx
                channelDataFallingAsleep = [channelDataFallingAsleep,EEG_table_timed{sec, eegChannels{ch}}{1,1}'];
                
            end

            if isempty(channelDataFallingAsleep)
               art_ch(ch,j) = ones(1,1);    
            else
               rms_all(ch, j) = rms(channelDataFallingAsleep);
            end 

        end
    else
        break; % Stop if the end index goes beyond the table
    end

end


%% Perform artefact detection based on selected method

for ch = 1:num_eeg_ch
    for ite = 1: iter % Number of iterations
        
        % If using RMS thresholding 
        if strcmp(method, 'RMS')
            med_rms = median(rms_all(ch, :),'omitnan');
            th = th_coeff*med_rms;
        
        % If using Otsu thresholding 
        elseif strcmp(method, 'Otsu')
        
            % Remove NaN values, then rescale to [0, 1]
        valid_rms = rms_all(~isnan(rms_all(ch, :)));
        scaled_rms = (valid_rms - min(valid_rms)) / (max(valid_rms) - min(valid_rms));
        otsu_level = graythresh(scaled_rms);
        th = otsu_level * (max(valid_rms) - min(valid_rms)) + min(valid_rms);
        else
            error('Unknown method. Use "RMS" or "Otsu".');
        end

        for i = 1:num_epc
            if rms_all(ch, i)> th
                num_epcs_rejected = num_epcs_rejected + 1;
                art_all(i) = 1;
                art_ch(ch,i) = 1;
                rms_all(ch, i) = NaN;
            end
        end
    end
end 


proportion_epochs_remaining = 1 - sum(art_all)/num_epc;
disp(['% Epochs remaining after rejecting ones ', num2str(th_coeff), ' times the thresholding method (', method, '): ', num2str(proportion_epochs_remaining)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
