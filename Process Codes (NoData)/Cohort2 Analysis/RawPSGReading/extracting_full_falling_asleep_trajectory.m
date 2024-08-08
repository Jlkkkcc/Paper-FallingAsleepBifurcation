
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Anastasia Ilina

% Preprocesses the patient .edf files and corresponding hypnograms. 
% Populates patient subtables with epoched EEG data for full sleep onset.
% Creates epoch-wise artefact masks for all different artefact thresholds
% (RMS method).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fixing directories 

% Add necessary directories

addpath('/home/anastasia/Dropbox/ISN work/Work on new dataset')
addpath('/mnt/LongTermStorage/Airforce_Surrey_psa/af0165')
addpath ('/home/anastasia/Dropbox/ISN work/Whole Code Pipelines/My code/edf-viewer-master')
addpath('/mnt/LongTermStorage/Airforce_Surrey_psa/HYPNOGRAMS')
addpath('/mnt/LongTermStorage/anastasia/Airforce_analysed')

% change directory to the folder with the hypnograms 
cd('/mnt/LongTermStorage/Airforce_Surrey_psa')

%% Get a list of all of the participants in the directory

% Get a list of all files and folders in the current directory
files = dir;

% Filter out all the files (keep only folders)
folders = files([files.isdir]);

% Remove the '.' and '..' directories
folders = folders(~ismember({folders.name}, {'.', '..', 'HYPNOGRAMS'}));

% Extract the folder names
participantNames = {folders.name};

% Define the path to the .mat file
filePath = '/mnt/LongTermStorage/anastasia/Airforce_analysed/allParticipantsTable.mat';

% Load the table
loadedData = load(filePath);

% Access the table
masterTable = loadedData.masterTable;

% How many minute to add to the falling asleep trace after Sleep Onset
Time2ExtractAfterOnset = 10; % In minutes
Time2ExtractAfterOnsetSec = Time2ExtractAfterOnset * 60;

th_coeff_list = [2, 2.5, 3, 3.5, 4];
iter = 2;

prop_epc_reject = NaN(length(participantNames), length(th_coeff_list));


if ~ismember('EpochedFallingAsleep', masterTable.Properties.VariableNames)
    masterTable.EpochedFallingAsleep = cell(height(masterTable), 1);
end

if ~ismember('ArtefactMask4Epochs', masterTable.Properties.VariableNames)
    masterTable.ArtefactMask4Epochs = cell(height(masterTable), 1);
end 

% Check and convert Night_Number to double
if iscell(masterTable.Night_Number)
    masterTable.Night_Number = cellfun(@str2double, masterTable.Night_Number);
elseif isstring(masterTable.Night_Number) 
    masterTable.Night_Number = str2double(masterTable.Night_Number);
end

% Ensure ifEmpty is logical
if iscell(masterTable.ifEmpty)
    masterTable.ifEmpty = cellfun(@(x) isequal(x, true), masterTable.ifEmpty);
end

% Define logical conditions for the entire masterTable
%isControlSleep = strcmp(masterTable.Type_of_Sleep_Night, 'Control');
isNormalBaseline = strcmp(masterTable.Baseline, 'Normal') | strcmp(masterTable.Baseline, 'Baseline');
isNightNumberInRange = masterTable.Night_Number >= 1 & masterTable.Night_Number <= 7;
isNotEmpty = ~masterTable.ifEmpty;

% Combine all conditions for the entire masterTable
validRows =  isNormalBaseline & isNightNumberInRange & isNotEmpty;

% Create a subset of masterTable based on conditions
subsetTable = masterTable(validRows, :);

% Process EEG data to create 6-second epochs with 3-second steps
eegSampleRate = 256; % 256 samples per second
epochDurationSec = 6; % 6 seconds
stepSec = 3; % 3 seconds

% Loop through the participant folders to extract the data from the hypnograms 
for i = 1:length(participantNames)
    participantName = participantNames{i};
    disp(participantName)

    % Open participant's folder 
    %ParticipantFolderDir = ['/mnt/LongTermStorage/Airforce_Surrey_psa/', participantName];
    %cd(ParticipantFolderDir);

    % Define logical condition for current patient within the subset
    %isCurrentPatient = strcmp(subsetTable.Patient_ID, participantName);

    % Define a participant-specific subset table
    %participantSubset = subsetTable(isCurrentPatient, :);


    % Construct the path to the .mat file
    savePathParticipantSubset = fullfile('/mnt/LongTermStorage/anastasia/Airforce_analysed/Participants_Epoched_EEG/', ['EpochedEEG_', participantName, '.mat']);
    
    % Check if the file exists before loading
    if exist(savePathParticipantSubset, 'file')
        % Load the table
        loadedData = load(savePathParticipantSubset);
        
        % Access the table (assuming 'participantSubset' is the variable name inside the MAT-file)
        participantSubset = loadedData.participantSubset;
        
        % Now you can use 'participantSubset' as needed
    else
        disp(['File not found: ', savePathParticipantSubset]);
    end
    

    % Iterate through the valid rows of the subset for current patient
    for n = 1:height(participantSubset)
        % Get the filename from the current row of the subset
        filename = fullfile('/mnt/LongTermStorage/Airforce_Surrey_psa/', [participantName, '/',  participantSubset.Filename{n}, '.edf']);
        disp(filename)
    
        % Check if the file exists
        if exist(filename, 'file')
            % Apply your function to read the EDF file
            [data, annotations] = edfread(filename);

            % Process EEG data to create epochs
            % Assuming data has columns: Time (in seconds), EEG channels
            eegChannels = data.Properties.VariableNames(1:end); % Excluding time column
            
        
            % Convert epoch indices to time in seconds
            bedtimeStartSec = participantSubset.AdjustedStartIndex(n) * 30;
            sleepOnsetStartSec = (participantSubset.SleepOnsetIndex(n) + participantSubset.AdjustedStartIndex(n))* 30 ;
            totalTimeSec = sleepOnsetStartSec + Time2ExtractAfterOnsetSec - bedtimeStartSec;
        
            % Extract EEG data
            FalingAsleepEEG = data(bedtimeStartSec: min(height(data), sleepOnsetStartSec+Time2ExtractAfterOnsetSec), :);
        
            % Calculate number of epochs
            numEpochs = floor((totalTimeSec - epochDurationSec) / stepSec) + 1;
            
            %% Detect epochs with artefacts
            % Create a table to store different coeff masks 
            % Create an empty table with column names 'Threshold', 'Artifact_Mask', and 'Percentage_Removed'
            artifact_table = table('Size', [0, 5], ...
                                    'VariableTypes', {'string', 'double', 'cell', 'cell', 'double'}, ...
                                    'VariableNames', {'Method', 'Threshold', 'Artifact_Mask_All', 'Artifact_Mask_Per_Channel',  'Percentage_Remaining_Pre_Sleep'});
            
            for th_idx = 1:length(th_coeff_list)
                
                th_coeff =th_coeff_list(th_idx);
                if th_coeff ~= 0
                    method = 'RMS';
                    [art_all, art_ch, proportion_epochs_remaining] = RMS_artefact_rejection(FalingAsleepEEG, numEpochs, epochDurationSec, stepSec, eegSampleRate, iter, th_coeff, method);
                    prop_epc_reject(i, th_idx) = proportion_epochs_remaining;
            
                  
                else
                    method = 'Otsu';
                    [art_all, art_ch, proportion_epochs_remaining] = RMS_artefact_rejection(FalingAsleepEEG, numEpochs, epochDurationSec, stepSec, eegSampleRate, iter, th_coeff, method);
                    prop_epc_reject(i, th_idx) = proportion_epochs_remaining;
                end 
                artifact_table = [artifact_table; {method, th_coeff, art_all, art_ch, proportion_epochs_remaining}];
            end 

            participantSubset.ArtifactMask4Epohcs{n} = artifact_table;


            %% Epoch the EEG data 

            % Initialize empty tables with all the necessary columns.
            % Variable `eegChannels` that contains the names of your EEG channels.
            allVariableNames = ['Timestamp', eegChannels];
            numVariables = length(allVariableNames);
            emptyCells = cell(0, numVariables);

            % Initialize the tables with the correct number of columns, filled with NaNs for numeric data or empty cells for non-numeric data
            % Initialize the tables with the correct number of columns
            eEpochsFallingAsleep = table('Size', [0, length(allVariableNames)], 'VariableTypes', ['double', repmat({'cell'}, 1, numVariables - 1)], 'VariableNames', allVariableNames);


            for j = 1:numEpochs
                % Calculate the starting index for this epoch in terms of rows in the 'data' table
                startTimeIdx = (j - 1) * stepSec + 1;
                % Calculate the ending index for this epoch in terms of rows in the 'data' table
                endTimeIdx = startTimeIdx + epochDurationSec - 1;

                % Debug: Print the current epoch number
                %disp(['Processing epoch: ', num2str(j)]);

                % Make sure the indices don't exceed the data limits
                if endTimeIdx <= height(FalingAsleepEEG)
                    % Initialize a matrix to hold the EEG data for this epoch
                    epochDataFallingAsleep = zeros(epochDurationSec * eegSampleRate, length(eegChannels));
                    % Extract the data for each channel for this epoch
                    for ch = 1:length(eegChannels)
                        channelDataFallingAsleep = [];

                        % Extract the EEG data from consecutive cells in the table and concatenate them
                        for sec = startTimeIdx:endTimeIdx
                            channelDataFallingAsleep = [channelDataFallingAsleep,FalingAsleepEEG{sec, eegChannels{ch}}{1,1}'];

                        end
                        % Place the concatenated channel data into the corresponding column of epochData
                        epochDataFallingAsleep(:, ch) = channelDataFallingAsleep;

                    end

                    % Add the epoch data to the epoch tables
                    eEpochsFallingAsleep(j, 'Timestamp') = {(j - 1) * stepSec};


                    for ch = 1:length(eegChannels)
                        % Debug: Print the current channel name
                        %disp(['Adding data to channel: ', eegChannels{ch}]);

                        eEpochsFallingAsleep{j, eegChannels{ch}} = {epochDataFallingAsleep(:, ch)};
                    end
                else
                    break; % Stop if the end index goes beyond the table
                end
            end


            % Store the processed epochs in the participant-specific subset table
            participantSubset.EpochedFallingAsleep{n} = eEpochsFallingAsleep;

            %% Plot rejected EEG epochs and save it in a separate folder
            night_number = join([participantSubset.Type_of_Sleep_Night(n), "_", participantSubset.Baseline(n), "_Night_", participantSubset.Night_Number(n)], '');
            saving_path = '/mnt/LongTermStorage/anastasia/Airforce_analysed/Participants_EEG_cleaning/';
            plot_rejected_epochs_falling_asleep(artifact_table, FalingAsleepEEG, participantName, 0, saving_path, night_number)
            

        else
            disp(['File not found: ', filename]);
            % Mark this row as empty
            %subsetTable.ifEmpty(n) = true;
            participantSubset.ifEmpty(n) = true;
            % Skip further processing for this row
            continue;
        end 
    end
   
    % Save the participant-specific subset table using MAT-file version 7.3
    savePathParticipantSubset = fullfile('/mnt/LongTermStorage/anastasia/Airforce_analysed/Participants_Epoched_EEG/Final_Participant_Files/', ['ArtifactFtEEGFallingAsleep_', participantName, '.mat']);

    % Use the '-v7.3' flag to save larger data
    save(savePathParticipantSubset, 'participantSubset', '-v7.3');
end
