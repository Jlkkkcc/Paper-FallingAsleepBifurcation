%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Anastasia Ilina

% Extracts the features from the pre-processed and epoched EEG data
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




% Assume participantSubset is your current table with EEG data
Fs = 256; % Sampling rate
python_directory = ''; % Your Python directory
iflinux = true; % Your operating system flag

% Initialize tables to store features for each epoch of each channel
numFeatures = 49;


SelectedChannels = {'F4_A1', 'C4_A1', 'P4_A1', 'O2_A1', 'LOC_A2', 'ROC_A1', 'F3_A2', 'C3_A2', 'P3_A2', 'O1_A2'};

% Loop through the participant folders to extract the data from the hypnograms 
for i = 1:length(participantNames)

    participantName = participantNames{i};
    disp(['Feature calculation started for Subject No. ',participantName])

    % Construct the path to the .mat file
    savePathParticipantSubset = fullfile('/mnt/LongTermStorage/anastasia/Airforce_analysed/Participants_Epoched_EEG/', ['ArtifactFtEEGFallingAsleep_', participantName, '.mat']);
    
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
    
    TableVariableNames =  participantSubset.EpochedFallingAsleep{1}.Properties.VariableNames;
   
    ChannelToExclude = 'EMG1_EMG2';
    excludeIndex = ismember(TableVariableNames, SelectedChannels);
    TableVariableNames = ['Timestamp', TableVariableNames(excludeIndex)];
   

    numColumns = length(TableVariableNames);

    %TableVariableTypes =  participantSubset.EpochedBedtimeBaseline5min{1}.Properties.VariableTypes;
    %TableVariableTypes = TableVariableTypes(~excludeIndex);

    
    % Loop through each epoch of each channel

    for n = 1:height(participantSubset)

        disp(['Calculating for night ',participantSubset.Filename{n}, '...'])
        disp('.... Falling Asleep Trajectory')
        
        if isempty(participantSubset.EpochedFallingAsleep{n})
            continue
        end
    
        participantSubset.FallingAsleepFeatures{n} = table('Size', [height(participantSubset.EpochedFallingAsleep{n}), numColumns], 'VariableTypes', ['double', repmat({'cell'}, 1, numColumns - 1)], 'VariableNames', TableVariableNames);
      


        timepoints = participantSubset.EpochedFallingAsleep{n}.Timestamp;
        participantSubset.FallingAsleepFeatures{n}.Timestamp = timepoints;
       
        % For bedtime baseline
        
        for j = 1:height(participantSubset.EpochedFallingAsleep{n})
            timepoint = participantSubset.EpochedFallingAsleep{n}(j,:);
            % Assuming your epochData is now a 2D array where columns are channels and rows are time points
            % Process each channel

            for ch = 2:(length(SelectedChannels)+1)
                currentChannel = TableVariableNames{ch};
                channelData = timepoint.(currentChannel){1,1};
                
                continuousZerosLength = 10;
                flatlineThreshold = 10;  % Number of consecutive points to check for flatlining


                % Check for continuous zeros (masked artefact)
                continuousZerosMask = conv(double(channelData == 0), ones(1, continuousZerosLength), 'same') == continuousZerosLength;
                
                % Check for flatlining
                flatlinedMask = conv(double(diff(channelData) == 0), ones(1, flatlineThreshold), 'same') >= (flatlineThreshold - 1);

            
                if any(continuousZerosMask) || any(flatlinedMask)
                    participantSubset.FallingAsleepFeatures{n}{j, ch} = {NaN(49,1)};
                else
                    [ft, ft_dscrp] = sleep_ft_extract_new(channelData, Fs, python_directory, iflinux);
                    % Store the features in the corresponding table
                    participantSubset.FallingAsleepFeatures{n}{j, ch} = {cell2mat(struct2cell(ft))};
                end
            end
        end
        
        % Save the participant-specific subset table using MAT-file version 7.3
        savePathParticipantSubset = fullfile('/mnt/LongTermStorage/anastasia/Airforce_analysed/Participants_Epoched_EEG/', ['ArtifactFtEEGFallingAsleep_', participantName, '.mat']);

        % Use the '-v7.3' flag to save larger data
        save(savePathParticipantSubset, 'participantSubset', '-v7.3');
    end
    
    disp(['Subject No. ',participantName ,' successfully calculated'])
   
end 