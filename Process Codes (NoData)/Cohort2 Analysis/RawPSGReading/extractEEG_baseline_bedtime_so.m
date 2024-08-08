
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Anastasia Ilina

% Preprocesses the patient .edf files and corresponding hypnograms. 
% Creates patient subtables where epoched EEG data for 5 minutes at the start
% of the bedtime and 5 minutes after sleep onset is stored.
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

% Add new columns to masterTable if they don't exist
%if ~ismember('BedtimeBaseline5min', masterTable.Properties.VariableNames)
%    masterTable.BedtimeBaseline5min = cell(height(masterTable), 1);
%end
%if ~ismember('SleepOnset5min', masterTable.Properties.VariableNames)
%    masterTable.SleepOnset5min = cell(height(masterTable), 1);
%end
if ~ismember('EpochedBedtimeBaseline5min', masterTable.Properties.VariableNames)
    masterTable.EpochedBedtimeBaseline5min = cell(height(masterTable), 1);
end
if ~ismember('EpochedSleepOnset5min', masterTable.Properties.VariableNames)
    masterTable.EpochedSleepOnset5min = cell(height(masterTable), 1);
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
    ParticipantFolderDir = ['/mnt/LongTermStorage/Airforce_Surrey_psa/', participantName];
    cd(ParticipantFolderDir);

    % Define logical condition for current patient within the subset
    isCurrentPatient = strcmp(subsetTable.Patient_ID, participantName);

    % Define a participant-specific subset table
    participantSubset = subsetTable(isCurrentPatient, :);

    % Iterate through the valid rows of the subset for current patient
    for n = 1:height(participantSubset)
        % Get the filename from the current row of the subset
        filename = [participantSubset.Filename{n}, '.edf'];
        disp(filename)
    
        % Check if the file exists
        if exist(filename, 'file')
            % Apply your function to read the EDF file
            [data, annotations] = edfread(filename);

            % Process EEG data to create epochs
            % Assuming data has columns: Time (in seconds), EEG channels
            eegChannels = data.Properties.VariableNames(1:end); % Excluding time column
            totalTimeSec = 5*60; % Total time in seconds of the EEG recording
        
            % Convert epoch indices to time in seconds
            bedtimeStartSec = participantSubset.AdjustedStartIndex(n) * 30;
            sleepOnsetStartSec = (participantSubset.SleepOnsetIndex(n) + participantSubset.AdjustedStartIndex(n)) * 30;
        
            % Extract EEG data
            durationSec = 5 * 60; % 5 minutes in seconds
            bedtimeBaseline5min = data(bedtimeStartSec + 1 : min(height(data), bedtimeStartSec + durationSec), :);
            sleepOnset5min = data(sleepOnsetStartSec + 1 : min(height(data), sleepOnsetStartSec + durationSec), :);
            
            %participantSubset.BedtimeBaseline5min{n} = bedtimeBaseline5min;
            %participantSubset.SleepOnset5min{n} = sleepOnset5min;
    
            % Initialize tables for epochs
            %eEpochsBedtime = table;
            %eEpochsSleepOnset = table;
        
            % Calculate number of epochs
            numEpochs = floor((totalTimeSec - epochDurationSec) / stepSec) + 1;
            

            % Initialize empty tables with all the necessary columns.
            % Variable `eegChannels` that contains the names of your EEG channels.
            allVariableNames = ['Timestamp', eegChannels];
            numVariables = length(allVariableNames);
            emptyCells = cell(0, numVariables);
            
            % Initialize the tables with the correct number of columns, filled with NaNs for numeric data or empty cells for non-numeric data
            % Initialize the tables with the correct number of columns
            eEpochsBedtime = table('Size', [0, length(allVariableNames)], 'VariableTypes', ['double', repmat({'cell'}, 1, numVariables - 1)], 'VariableNames', allVariableNames);
            eEpochsSleepOnset = table('Size', [0, length(allVariableNames)], 'VariableTypes', ['double', repmat({'cell'}, 1, numVariables - 1)], 'VariableNames', allVariableNames);


            for j = 1:numEpochs
                % Calculate the starting index for this epoch in terms of rows in the 'data' table
                startTimeIdx = (j - 1) * stepSec + 1;
                % Calculate the ending index for this epoch in terms of rows in the 'data' table
                endTimeIdx = startTimeIdx + epochDurationSec - 1;
                
                % Debug: Print the current epoch number
                %disp(['Processing epoch: ', num2str(j)]);
            
                % Make sure the indices don't exceed the data limits
                if endTimeIdx <= height(bedtimeBaseline5min)
                    % Initialize a matrix to hold the EEG data for this epoch
                    epochDataBaseline = zeros(epochDurationSec * eegSampleRate, length(eegChannels));
                    epochDataSleepOnset = zeros(epochDurationSec * eegSampleRate, length(eegChannels));
                    % Extract the data for each channel for this epoch
                    for ch = 1:length(eegChannels)
                        channelDataBaseline = [];
                        channelDataSleepOnset = [];

                        % Extract the EEG data from consecutive cells in the table and concatenate them
                        for sec = startTimeIdx:endTimeIdx
                            channelDataBaseline = [channelDataBaseline, bedtimeBaseline5min{sec, eegChannels{ch}}{1,1}'];
                            channelDataSleepOnset = [channelDataSleepOnset, sleepOnset5min{sec, eegChannels{ch}}{1,1}'];
                        end
                        % Place the concatenated channel data into the corresponding column of epochData
                        epochDataBaseline(:, ch) = channelDataBaseline;
                        epochDataSleepOnset(:, ch) = channelDataSleepOnset;
                    end
            
                    % Add the epoch data to the epoch tables
                    eEpochsBedtime(j, 'Timestamp') = {(j - 1) * stepSec};
                    eEpochsSleepOnset(j, 'Timestamp') = {(j - 1) * stepSec};
            
                    for ch = 1:length(eegChannels)
                        % Debug: Print the current channel name
                        %disp(['Adding data to channel: ', eegChannels{ch}]);
            
                        eEpochsBedtime{j, eegChannels{ch}} = {epochDataBaseline(:, ch)};
                        eEpochsSleepOnset{j, eegChannels{ch}} = {epochDataSleepOnset(:, ch)};
                    end
                else
                    break; % Stop if the end index goes beyond the table
                end
            end


            % Store the processed epochs in the participant-specific subset table
            participantSubset.EpochedBedtimeBaseline5min{n} = eEpochsBedtime;
            participantSubset.EpochedSleepOnset5min{n} = eEpochsSleepOnset;
   
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
    savePathParticipantSubset = fullfile('/mnt/LongTermStorage/anastasia/Airforce_analysed/Participants_Epoched_EEG/', ['EpochedEEG_', participantName, '.mat']);

    % Use the '-v7.3' flag to save larger data
    save(savePathParticipantSubset, 'participantSubset', '-v7.3');
end

% Define the path for saving the subset
%savePathSubset = '/mnt/LongTermStorage/anastasia/Airforce_analysed/allParticipantsTableSubsetWithEEG.mat';

% Save the subsetTable
%save(savePathSubset, 'subsetTable');

% Save it as .csv as well in case we want to use python later on.
%writetable(masterTable, strrep(savePath, '.mat', '.csv'));
