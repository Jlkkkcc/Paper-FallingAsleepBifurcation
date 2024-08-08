%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Anastasia Ilina

% Uses the artefact mask to swap the features in the epochs with artefacts
% with NaNs. Masks with different artefact thresholds can be applied 
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

% How many minute to add to the falling asleep trace after Sleep Onset
Time2ExtractAfterOnset = 10; % In minutes
Time2ExtractAfterOnsetSec = Time2ExtractAfterOnset * 60;

ChosenArtifactThreshold = 3;



% Loop through the participant folders to extract the data from the hypnograms 
for i = 1:length(participantNames)
    participantName = participantNames{i};
    disp(participantName)

     % Save the participant-specific subset table using MAT-file version 7.3
    savePathParticipantSubset = fullfile('/mnt/LongTermStorage/anastasia/Airforce_analysed/Participants_Epoched_EEG/Final_Participant_Files/', ['ArtifactFtEEGFallingAsleep_', participantName, '.mat']);
    
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
        
    if ~ismember('CleanFallingAsleepFeatures', participantSubset.Properties.VariableNames)
        participantSubset.CleanFallingAsleepFeatures = cell(height(participantSubset), 1);
    end

    % Iterate through the valid rows of the subset for current patient
    for n = 1:height(participantSubset)
        disp(['Cleaning fts for night ',participantSubset.Filename{n}, '...'])
        
        % check if the artifact table for the night is present
        if isempty(participantSubset.ArtifactMask4Epohcs{n})
           continue
        end 
        
        % Extract the artifact table for the night
        artifact_table = participantSubset.ArtifactMask4Epohcs{n};

        % Select an artifact mask according to the chosen threshold
        selectedArtifactMask = artifact_table.Artifact_Mask_Per_Channel(artifact_table.Threshold == ChosenArtifactThreshold);
        selectedArtifactMask = selectedArtifactMask{1};
        
        % Initialise clean feature table with raw feature values
        cleanFts = participantSubset.FallingAsleepFeatures{n};
        
        % Go through each channel (column), and set the feature values of
        % epochs with artifact to NaN
        
        for col = 1:(size(cleanFts, 2) - 1)
            if strcmp(cleanFts.Properties.VariableNames{1}, 'Timestamp')
                clean_ft_ch  = cleanFts{:, col+1};
            else
                clean_ft_ch  = cleanFts{:, col};
            end 
            
            % extract artifact mask for the current channel
            artifact_ch = squeeze(selectedArtifactMask(col, :));
            
            % check if there are any artifacts for this channel
            if sum(artifact_ch) > 0
                clean_ft_ch(artifact_ch == 1) = {NaN(49, 1)};
            end
            
            if strcmp(cleanFts.Properties.VariableNames{1},'Timestamp')
                cleanFts{:, col+1} = clean_ft_ch;
            else
                cleanFts{:, col} = clean_ft_ch;
            end
        end 
   
            
        % Save it to the participant table for this night 
        participantSubset.CleanFallingAsleepFeatures{n} = cleanFts;

    end
   
    % Save the participant-specific subset table using MAT-file version 7.3
    savePathParticipantSubset = fullfile('/mnt/LongTermStorage/anastasia/Airforce_analysed/Participants_Epoched_EEG/Final_Participant_Files/', ['ArtifactFtEEGFallingAsleep_', participantName, '.mat']);

    % Use the '-v7.3' flag to save larger data
    save(savePathParticipantSubset, 'participantSubset', '-v7.3');
end
