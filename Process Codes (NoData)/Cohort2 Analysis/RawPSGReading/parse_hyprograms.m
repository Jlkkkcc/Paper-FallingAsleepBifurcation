
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Anastasia Ilina

% Parses the hypnogram text files of the patien to create a table containing a summary
% all the sleep information (hypnograms, information about a given night of
% sleep, etc). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fixing directories 

% Add necessary directories

addpath('/home/anastasia/Dropbox/ISN work/Work on new dataset')
addpath('/mnt/LongTermStorage/Airforce_Surrey_psa/af0165')
addpath ('/home/anastasia/Dropbox/ISN work/Whole Code Pipelines/My code/edf-viewer-master')
addpath('/mnt/LongTermStorage/Airforce_Surrey_psa/HYPNOGRAMS')
addpath('/mnt/LongTermStorage/anastasia/Airforce_analysed')

% change directory to the folder with the hypnograms 
cd('/mnt/LongTermStorage/Airforce_Surrey_psa/HYPNOGRAMS')

%% Get a list of all of the participants in the directory

% Get a list of all files and folders in the current directory
files = dir;

% Filter out all the files (keep only folders)
folders = files([files.isdir]);

% Remove the '.' and '..' directories
folders = folders(~ismember({folders.name}, {'.', '..'}));

% Extract the folder names
participantNames = {folders.name};

% Display the folder names
%disp(folderNames);

% Initialize the master table with an additional column for Filename
columnNames = {'Filename', 'Patient_ID', 'Type_of_Sleep_Night', 'Baseline', 'Night_Number', 'Recording_Date', 'Start_of_Recording', 'Recording_Duration', 'SleepStagesInt', 'SleepStagesString'};
variableTypes = {'string', 'string', 'string', 'string', 'string', 'string', 'string', 'int32', 'cell', 'cell'};
masterTable = table('Size', [0, length(columnNames)], 'VariableTypes', variableTypes, 'VariableNames', columnNames);

% Loop through the participant folders to extract the data from the
% hypnograms 
for i = 1:length(participantNames)
    participantName = participantNames(i);
    participantName = participantName{1,1};
    disp(participantName)
    
    % Open participant's folder 
    ParticipantFolderDir = ['/mnt/LongTermStorage/Airforce_Surrey_psa/HYPNOGRAMS/', participantName];
    cd(ParticipantFolderDir);

    % Get a list of all .txt files in the current directory
    txtFiles = dir('*.txt');

    % Initialize a table for the participant
    %columnNames = {'Patient_ID', 'Type_of_Sleep_Night', 'Baseline', 'Night_Number', 'Recording_Date', 'Start_of_Recording', 'Recording_Duration', 'SleepStagesInt', 'SleepStagesString'};
    %variableTypes = {'string', 'string', 'string', 'string', 'string', 'string', 'int32', 'cell', 'cell'}; 

    % Initialize a table for the participant with an additional column for Filename
    participantTable = table('Size', [0, length(columnNames)], 'VariableTypes', variableTypes, 'VariableNames', columnNames);

    
    % Iterate through each file
    for i = 1:length(txtFiles)
        filename = txtFiles(i).name;
        
        
        % Open the file
        fileID = fopen(filename, 'r');
        
        % Check if the file opened successfully
        if fileID == -1
            warning('Unable to open file: %s', filename);
            continue; % Skip to the next file
        end
    
        % Parse the filename
        tokens = regexp(filename, 'AFOSR_AF(\d+)_([S][ER][BER][NR]\d+)_(\d\d[A-Z]{3}\d\d)_(\d\d\d\d)\.txt', 'tokens');
        if isempty(tokens)
            continue; % Skip if the filename doesn't match the pattern
        end
        tokens = tokens{1};
        patientID = ['af', lower(tokens{1})];
        sleepTypeLetter = tokens{2}(2); % 'E' for normal, 'R' for restricted
        
        % Determine sleep type
        if sleepTypeLetter == 'E'
            sleepType = "Control";
        elseif sleepTypeLetter == 'R'
            sleepType = "Restricted";
        else
            sleepType = "Unknown"; % Fallback case
        end
        
        % Determine whether it's a baseline recording
        baselineTypeLetter = upper(tokens{2}(3));
        if baselineTypeLetter == 'B'
            baseline = 'Baseline';
        else
            baseline = 'Normal';
        end 

        % Determine the rest

        nightNumber = tokens{2}(end);
        if nightNumber == '8'
            baseline = 'Recovery';
        end

        recordingDate = tokens{3};
        startTime = insertAfter(tokens{4}, 2, ':');
        
        % Initialize arrays for storing sleep stages
        sleepStages = []; % For integer representation
        sleepStagesCell = {}; % For string representation

        % Construct the filename to be saved in the table
        % It includes everything up to and including the recording date
        filenameToSave = ['AFOSR_AF', tokens{1}, '_', tokens{2}];
    
        % Read file line by line
        while ~feof(fileID)
            line = fgetl(fileID); % Read one line from the file
            
            if ischar(line) % Check if the line is not the end of the file
                % Convert symbol to integer representation
                switch line
                    case '?'
                        stageInt = NaN;
                        stageStr = '';
                    case 'W'
                        stageInt = 0;
                        stageStr = 'W';
                    case 'R'
                        stageInt = 5;
                        stageStr = 'R';
                    otherwise
                        if any(line == '1234') % Check if line is '1', '2', '3', or '4'
                            stageInt = str2double(line);
                            stageStr = ['N', line]; % N1, N2, N3, N4
                        else
                            continue; % Skip lines with unexpected symbols
                        end
                end
    
                % Append to the arrays
                sleepStages(end+1) = stageInt;
                sleepStagesCell{end+1} = stageStr;
            end
        end
        
        % Calculate the recording duration
        recordingDuration = length(sleepStages) / 2;

        % Add the data to the table
         % Add the data to the table
        newRow = {filenameToSave, patientID, sleepType, baseline, nightNumber, recordingDate, startTime, recordingDuration, {sleepStages}, {sleepStagesCell}};
        participantTable = [participantTable; newRow];
    
        % Close the file
        fclose(fileID);

    end
    masterTable = [masterTable; participantTable];
    % Display the table
    %disp(participantTable);

end

% Define the path
savePath = '/mnt/LongTermStorage/anastasia/Airforce_analysed/allParticipantsTable.mat';

% Save the table
save(savePath, 'masterTable');

% Save it as .csv as well in case we want to use python later on.
writetable(masterTable, strrep(savePath, '.mat', '.csv'));

