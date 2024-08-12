function plot_rejected_epochs_falling_asleep(artifact_table, falling_asleep_eeg, sbj_id, ifallchmask, saving_path, night_number)

% Replace 'eeg_data' and 'artifact_masks' with your actual variables containing the EEG data and artifact masks, respectively.
% Make sure the dimensions of 'eeg_data' and 'artifact_masks' are compatible.


% Sampling rate of the EEG data (replace with your actual sampling rate)
sampling_rate = 256;  % In Hz
epoch_duration = 6; % In seconds
overlap = 3;
    
% Length of the EEG recording in seconds
recording_duration = height(falling_asleep_eeg);
% Time vector corresponding to the EEG data
time_axis = linspace(0, recording_duration, recording_duration*sampling_rate);

eegChannels = falling_asleep_eeg.Properties.VariableNames;
num_chs = length(falling_asleep_eeg.Properties.VariableNames);
thresholds = artifact_table.Threshold;



% Assuming sbj_id is a variable containing the subject ID
folder_name = join([sbj_id, '/', night_number, '/'],'');
if ifallchmask
    subfolder_name = 'intersection_of_channels_filter';
else
    subfolder_name = 'filter_per_channel';
end 
saving_path_sbj = join([saving_path, folder_name], '');

% Check if the folder exists
if ~exist(saving_path_sbj, 'dir')
    % If the folder does not exist, create it
    mkdir(saving_path_sbj);
    disp(['Folder "', saving_path_sbj, '" created.']);
end

% Change to the directory of the folder
%cd(folder_name);
%disp(['Current directory changed to "', folder_name, '".']);

saving_path_subfolder = join([saving_path_sbj,'/', subfolder_name], '');
% Check if the folder exists
if ~exist(saving_path_subfolder, 'dir')
    % If the folder does not exist, create it
    mkdir(saving_path_subfolder);
    disp(['Folder "', saving_path_subfolder, '" created.']);
end

% Change to the directory of the folder
%cd(subfolder_name);
%disp(['Current directory changed to "', subfolder_name, '".']);


eeg_mat = NaN(num_chs, recording_duration*sampling_rate);
for ch = 1:num_chs
    channelDataFallingAsleep = [];
    
    % Extract the EEG data from consecutive cells in the table and concatenate them
    for sec = 1:recording_duration
        channelDataFallingAsleep = [channelDataFallingAsleep,falling_asleep_eeg{sec, eegChannels{ch}}{1,1}'];
       
    end
    eeg_mat(ch, :) = channelDataFallingAsleep;

end 
    
for ch = 1:num_chs
    eeg_mat_ch = eeg_mat(ch, :);
    
    
    for th = 1:length(thresholds)
        cur_th = thresholds(th);
        
        if ifallchmask
            art_mask = artifact_table(th, :).Artifact_Mask_All{1,1};
            

        else
            art_mask = artifact_table(th, :).Artifact_Mask_Per_Channel{1,1}(ch,:);
           
        end
        num_epochs = max(size(art_mask));
        time_axis_epochs = linspace(0, num_epochs * (epoch_duration-overlap) - overlap, num_epochs);
        

        figure;
        
        % Plot the EEG data 
        subplot(1, 1, 1);
        plot(time_axis, eeg_mat_ch);
        xlabel('Time (s)');
        ylabel('Amplitude (µV)');
        title(['Channel: ', strrep(eegChannels{ch}, '_', '-'), '. Threshold: below ', num2str(cur_th), ' RMS']);
        %ylim([-0.4, 0.4]);
        xlim([0, recording_duration])

        % Highlight the rejected epochs in Pre-Sleep (set the color of the corresponding line segments to red)
        rejected_epochs = find(art_mask == 1);
        if ~isempty(rejected_epochs)
            hold on;
            % For Pre-Sleep
            for i = 1:length(rejected_epochs)
                if mod(i,2) == 0
                    continue
                else
                    start_time = time_axis_epochs(rejected_epochs(i));
                    end_time = start_time + epoch_duration;
                    % Find the corresponding time points within the rejected epoch
                    time_points_within_rejected_epoch = (time_axis > start_time) & (time_axis <= end_time);
                    
                    % Set the color of the line segments within the rejected epoch to red
                    plot(time_axis(time_points_within_rejected_epoch), eeg_mat_ch(time_points_within_rejected_epoch), 'r');
            
                end
            end 
            hold off;
        end

        figure_filename = join([saving_path_subfolder, '/', 'Channel_', num2str(ch), '_Threshold_', num2str(cur_th), '.jpg'], '');
        saveas(gcf, figure_filename, 'jpg');
        disp(join(['Figure saved as ', figure_filename],''));
        close(gcf);
    end
end