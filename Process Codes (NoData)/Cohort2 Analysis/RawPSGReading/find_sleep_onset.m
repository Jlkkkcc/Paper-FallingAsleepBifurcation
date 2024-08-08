%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Anastasia Ilina

% Processes patients' hypnogram data to clean it and identify sleep onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add necessary directories

addpath('/home/anastasia/Dropbox/ISN work/Work on new dataset')
addpath('/home/anastasia/Dropbox/ISN work/Work on new dataset/helper_functions')
addpath('/mnt/LongTermStorage/Airforce_Surrey_psa/af0165')
addpath ('/home/anastasia/Dropbox/ISN work/Whole Code Pipelines/My code/edf-viewer-master')
addpath('/mnt/LongTermStorage/Airforce_Surrey_psa/HYPNOGRAMS')
addpath('/mnt/LongTermStorage/anastasia/Airforce_analysed')

% change directory to the folder with the hypnograms 
cd('/mnt/LongTermStorage/anastasia/Airforce_analysed')

% Define the Sleep Onset period extraction parameters
timeMinutesBeforeSleepOnset = 30;
timeMinutesAfterSleepOnset = 5;
hypnogramSecondsEpochDuration = 30;

% Define the path to the .mat file
filePath = '/mnt/LongTermStorage/anastasia/Airforce_analysed/allParticipantsTable.mat';

% Load the table
loadedData = load(filePath);

% Access the table
masterTable = loadedData.masterTable;

% Number of rows in the table
numRows = height(masterTable);

% Define dynamic column names
postOnsetColumnName = sprintf('AfterOnset_%d_min', timeMinutesAfterSleepOnset);
preOnsetColumnName = sprintf('BeforeOnset_%d_min', timeMinutesBeforeSleepOnset);

% Add new columns to masterTable for sleep onset-related information
% Add new columns to masterTable if they don't exist
if ~ismember('ifEmpty', masterTable.Properties.VariableNames)
    masterTable.ifEmpty = false(height(masterTable), 1);
end
if ~ismember('ifRecordingStartIsShifted', masterTable.Properties.VariableNames)
    masterTable.ifRecordingStartIsShifted = false(height(masterTable), 1);
end
if ~ismember('AdjustedStartIndex', masterTable.Properties.VariableNames)
    masterTable.AdjustedStartIndex = zeros(height(masterTable), 1);
end
if ~ismember('CleanSleepStagesInt', masterTable.Properties.VariableNames)
    masterTable.CleanSleepStagesInt = cell(height(masterTable), 1);
end
if ~ismember('SleepOnsetIndex', masterTable.Properties.VariableNames)
    masterTable.SleepOnsetIndex = NaN(height(masterTable), 1);
end
if ~ismember(postOnsetColumnName, masterTable.Properties.VariableNames)
    masterTable.(postOnsetColumnName) = cell(height(masterTable), 1);
end
if ~ismember(preOnsetColumnName, masterTable.Properties.VariableNames)
    masterTable.(preOnsetColumnName) = cell(height(masterTable), 1);
end
if ~ismember('ifCleanOnset', masterTable.Properties.VariableNames)
    masterTable.ifCleanOnset = false(height(masterTable), 1);
end
if ~ismember('Time2Sleep', masterTable.Properties.VariableNames)
    masterTable.Time2Sleep = NaN(height(masterTable), 1);
end



% Iterate through each row
for i = 1:numRows
    % Access data from the ith row
    currentRow = masterTable(i, :);
    
    patientID = currentRow.Patient_ID;
    sleepstagearray= currentRow.SleepStagesInt;
    sleepstagearray = sleepstagearray{1,1};

    if currentRow.Recording_Duration == 0
        masterTable.ifEmpty(i) = true;
        continue
    end

    [postOnsetScores, preOnsetScores, ifcleanonset, time2sleep, onsetIndex, adjustedStartIndex, ifRecordingStartIsShifted, sleepstagearrayClean] = sleepOnsetAnalysis(sleepstagearray, timeMinutesBeforeSleepOnset, timeMinutesAfterSleepOnset, hypnogramSecondsEpochDuration);
    
    masterTable.ifRecordingStartIsShifted(i) = ifRecordingStartIsShifted;
    masterTable.AdjustedStartIndex(i) = adjustedStartIndex;
    masterTable.CleanSleepStagesInt{i} = sleepstagearrayClean;
    masterTable.SleepOnsetIndex(i) = onsetIndex;
    masterTable.(postOnsetColumnName)(i) = {postOnsetScores};
    masterTable.(preOnsetColumnName)(i) = {preOnsetScores};
    masterTable.ifCleanOnset(i) = ifcleanonset;
    masterTable.Time2Sleep(i) = time2sleep;


end

% Save the updated masterTable
save(filePath, 'masterTable');

% Define the path for the .csv file
csvFilePath = '/mnt/LongTermStorage/anastasia/Airforce_analysed/allParticipantsTable_updated.csv';

% Save the updated masterTable as a .csv file
writetable(masterTable, csvFilePath);


