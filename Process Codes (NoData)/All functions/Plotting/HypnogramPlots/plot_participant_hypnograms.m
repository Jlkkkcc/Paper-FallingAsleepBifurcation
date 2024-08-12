%% Fixing directories 

% Add necessary directories

addpath('/home/anastasia/Dropbox/ISN work/Work on new dataset')
addpath('/mnt/LongTermStorage/Airforce_Surrey_psa/af0165')
addpath ('/home/anastasia/Dropbox/ISN work/Whole Code Pipelines/My code/edf-viewer-master')
addpath('/mnt/LongTermStorage/Airforce_Surrey_psa/HYPNOGRAMS')
addpath('/mnt/LongTermStorage/anastasia/Airforce_analysed')
addpath('/home/anastasia/Dropbox/ISN work/Whole Code Pipelines/My code/Shaded Error Bar')
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


SelectedChannels = {'F4_A1', 'C4_A1', 'P4_A1', 'O2_A1', 'LOC_A2', 'ROC_A1', 'F3_A2', 'C3_A2', 'P3_A2', 'O1_A2'};

sleepOnsetDuration = 60; %in seconds 

saveDirImg = '/home/anastasia/Dropbox/ISN work/Work on new dataset/Full Sleep Onset/Hypnograms';

epoch_dur = 0.5; % in minutes
timeAfterSleepOnset = 10; % in minutes

% Load the fitting_results_table

load('/mnt/LongTermStorage/anastasia/Airforce_analysed/Individual Bifurcation Fitting/New_Fitting_Results_30min.mat')


% Get rid of the duplicates 
% Sort the table by 'participantName' and 'nightNumber' in descending order to ensure the last occurrence comes first
sortedTable = sortrows(fitting_result_table, {'participantName', 'nightNumber'}, 'descend');

% Then find unique combinations of 'participantName' and 'nightNumber', keeping the first occurrence (which is actually the last occurrence in the original table)
[~, idx] = unique(sortedTable(:, {'participantName', 'nightNumber'}), 'rows', 'stable');

% Extract these rows to get a table without duplicates, preserving the last occurrence of each duplicate
fitting_result_table = sortedTable(idx, :);


% Add columns for sleep stats if they don't already exist
requiredColumns = {'EarliestNREM2PreSleep', 'LatestNREM2PreSleep', 'NumberNREM2PreSleep', 'TotalNREM1DurationPreSleep', 'TotalDurationPreSleep', 'ProportionNREM1PreSleep', 'AverageSleepDepthPreSleep','InconsistencyIndex', 'FluctuationIndex', 'RelativeWakeDurationPostSleep', 'RelativeNREM1DurationPostSleep', 'EarliestWakeTransitionPostSleep', 'AverageSleepDepthPostSleep', 'FluctuationIndexPostSleep', 'STD_S', 'Smoothness_S', 'Kurtosis_S', 'RMSSD_S'};
for colName = requiredColumns
    if ~any(strcmp(fitting_result_table.Properties.VariableNames, colName))
        fitting_result_table.(colName{1}) = NaN(height(fitting_result_table), 1); % Initialize with NaNs
    end
end

% Loop through the participant folders to extract the data from the hypnograms 
for i = 1:length(participantNames)

    participantName = participantNames{i};
   
    saveDirHypnogramParticipant = [saveDirImg, '/', participantName];

    if ~exist(saveDirHypnogramParticipant, 'dir')
        mkdir(saveDirHypnogramParticipant)
    end

    disp(['Hypnogram visualisation started for Subject No. ',participantName])
    
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

    for n = 1:8%height(participantSubset)
        fprintf('Night Number: %g\n', n);

        if isempty(participantSubset.CleanSleepStagesInt{n})
            continue
        end

        allNightHypnogram = participantSubset.CleanSleepStagesInt{n};
        sleepOnsetIndex = participantSubset.SleepOnsetIndex(n);
        S_variable = participantSubset.S_variable{n};

        numEpochsAfterSO = timeAfterSleepOnset/epoch_dur;
        
        if length(allNightHypnogram) >= sleepOnsetIndex+numEpochsAfterSO
            selectedHypnogram = allNightHypnogram(1:sleepOnsetIndex+numEpochsAfterSO);
            preSleepPeriod = allNightHypnogram(1:sleepOnsetIndex-1);
            postSleepPeriod = allNightHypnogram(sleepOnsetIndex + 1:sleepOnsetIndex+numEpochsAfterSO);   
        else
            selectedHypnogram = allNightHypnogram(1:end);
            preSleepPeriod = allNightHypnogram(1:sleepOnsetIndex-1);
            postSleepPeriod = allNightHypnogram(sleepOnsetIndex + 1:end);
        end 

        % %% Plot it 
        % timelineHypnogram = linspace((-sleepOnsetIndex+1)*epoch_dur, timeAfterSleepOnset, length(selectedHypnogram)); 
        % 
        % allNightHypnogram = participantSubset.CleanSleepStagesInt{n};
        % 
        % % Define sleep stage names
        % stageNames = {'Awake', 'NREM Stage 1', 'NREM Stage 2', 'NREM Stage 3', 'NREM Stage 4', 'REM'};
        % 
        % % Plot the hypnogram
        % 
        % night_fig = figure;
        % plot(timelineHypnogram, selectedHypnogram, 'b.-', 'LineWidth', 2);
        % 
        % % Set y-axis labels
        % set(gca, 'YTick', 0:5, 'YTickLabel', stageNames);
        % 
        % set(gca,'FontSize', 14)
        % set(gca,'TickDir','out')
        % set(gca,'ticklength',2*get(gca,'ticklength'))
        % set(gca,'lineWidth',2)
        % grid on;
        % 
        % % Set x and y-axis labels
        % xlabel('Time (minutes)');
        % %ylabel('Sleep Stage');
        % 
        % 
        % % Construct the title with line breaks
        % graphTitle = {'Sleep Stages Hypnogram', ...    
        %               ['Participant ', participantName, ' on Night ', num2str(n)],...
        %               ''};
        % 
        % % Set the title
        % title(graphTitle, 'Interpreter', 'none');
        % 
        % ax = gca;
        % 
        % ylim([-0.5,4.5])
        % 
        % % Draw a vertical line at time 0
        % line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r');
        % 
        % 
        % 
        % % Add a legend
        % %legend('Sleep Stages', 'Sleep Onset');
        % 
        % % Construct filenames for saving
        % jpegFilename = fullfile(saveDirHypnogramParticipant, ['Participant_', participantName, '_Night', num2str(n), '_Plot.jpg']);
        % figFilename = fullfile(saveDirHypnogramParticipant, ['Participant_', participantName, '_Night', num2str(n), '_Plot.fig']);
        % 
        % % Save the plot as JPEG
        % saveas(night_fig, jpegFilename);
        % 
        % % Save the figure as .fig file
        % savefig(night_fig, figFilename);
        % 
        % close(night_fig)
        
        %% Calculate Sleep Stats 
        

        % Assuming preSleepPeriod and SleepOnsetIndex are defined
        % epoch_dur is assumed to be defined, representing the duration of each epoch (e.g., 30 seconds)
        
        % Initialization
        ifPreSleepNREM2 = false;
        numberPreSleepNREM2 = 0;
        earliestNREM2BeforeSleep = NaN; % Initialize to NaN to handle cases with no NREM2
        latestNREM2BeforeSleep = NaN; % Initialize to NaN to handle cases with no NREM2
        
        if sum(preSleepPeriod == 2) > 0
            ifPreSleepNREM2 = true;
            numberPreSleepNREM2 = sum(preSleepPeriod == 2);
        
            % Correct use of find to get indices of NREM2
            nrem2Indices = find(preSleepPeriod == 2);
            
            % Correct calculation of distances
            earliestNREM2Index = nrem2Indices(1); % First occurrence of NREM2
            latestNREM2Index = nrem2Indices(end); % Last occurrence of NREM2
            earliestNREM2BeforeSleep = epoch_dur * (sleepOnsetIndex - earliestNREM2Index);
            latestNREM2BeforeSleep = epoch_dur * (sleepOnsetIndex - latestNREM2Index);
        end
        
        % Ensure you handle NaN cases or no NREM2 cases gracefully when displaying results
        fprintf('Earliest Failed Sleep Attempt (NREM2): %g minutes\n', earliestNREM2BeforeSleep);
        fprintf('Latest Failed Sleep Attempt (NREM2): %g minutes\n', latestNREM2BeforeSleep);
        fprintf('Number of Failed Sleep Attempt (NREM2): %d\n', numberPreSleepNREM2);


        % Example input
        %preSleepPeriod = [0 0 0 1 1 1 0 1 0 0 1 1 0 0];
        
        % Find NREM1 bouts starts
        startNREM1 = find(diff([0 preSleepPeriod]) == 1);
        
        % Initialize endNREM1 array
        endNREM1 = zeros(size(startNREM1));
        
        % For each NREM1 start, find the next change in sleep stage
        for ind = 1:length(startNREM1)
            nextChangeIndex = find(preSleepPeriod(startNREM1(ind):end) ~= 1, 1, 'first');
            if isempty(nextChangeIndex)
                % If no change is found, the end is the last index of preSleepPeriod
                endNREM1(ind) = length(preSleepPeriod);
            else
                % Adjust for the actual position
                endNREM1(ind) = startNREM1(ind) + nextChangeIndex - 2;
            end
        end
        
        % Calculate duration of each NREM1 bout and total NREM1 duration
        NREM1Durations = (endNREM1 - startNREM1 + 1) * 0.5; % Duration in minutes (each index is 30 seconds)
        totalNREM1Duration = sum(NREM1Durations); % Total duration in NREM1

        % Calculate total duration of preSleepPeriod
        totalDuration = length(preSleepPeriod) * 0.5; % Total duration in minutes
        
        % Calculate the proportion of time spent in NREM1
        proportionNREM1 = totalNREM1Duration / totalDuration;
        
        % Calculate inconsistency metric (simple version)
        % Count of transitions can indicate inconsistency
        numTransitions = sum(abs(diff(preSleepPeriod)) == 1); % Counting transitions
        inconsistencyIndex = numTransitions / totalNREM1Duration; % Adjusted by total NREM1 duration
        
        % Fluctuation Index
        % Simple measure based on the total number of state changes
        fluctuationIndex = numTransitions / totalDuration; % Adjusted by total duration

        % Compute sleep depth pre-onset
        preonsetSleepDepth = mean(preSleepPeriod);
        
        % Output the results
        fprintf('NREM1 Durations (minutes): '); fprintf('%g ', NREM1Durations); fprintf('\n');
        fprintf('Total Duration in NREM1: %g minutes\n', totalNREM1Duration);
        fprintf('Total Pre-sleep Duration: %g minutes\n', totalDuration);
        fprintf('Proportion of Time in NREM1: %g\n', proportionNREM1);
        fprintf('Inconsistency Index: %g\n', inconsistencyIndex);
        fprintf('Fluctuation Index: %g\n', fluctuationIndex);
        fprintf('Pre-Sleep-Onset Sleep Depth: %g\n', preonsetSleepDepth);


        %% Calculate Post-Falling Asleep Stats
        %postSleepPeriod = [2, 2, 3, 3, 2, 0, 2, 3, 0, 1, 0, 2, 5, 5, 0, 2, 2, 0, 3, 3];

        
        % Count wake (0) epochs in postSleepPeriod for relative duration calculations
        numWakeEpochs = sum(postSleepPeriod == 0);
        numNREM1Epochs = sum(postSleepPeriod == 1);
        
        % Total duration of wake epochs in minutes
        totalWakeDuration = numWakeEpochs * epoch_dur;
        totalNREM1DurationPost = numNREM1Epochs * epoch_dur;
        
        % Total duration of the postSleepPeriod
        totalPostDuration = length(postSleepPeriod) * epoch_dur;
        
        % Relative duration of wake epochs
        relativeWakeDuration = totalWakeDuration / totalPostDuration;
        relativeNREM1Duration = totalNREM1DurationPost / totalPostDuration;
        
        % Count of total transitions
        totalTransitions = sum(abs(diff(postSleepPeriod)) > 0);
        
        % Calculate the transitions to wake. A transition to wake is identified by a change from any non-zero state to 0.
        % Initialize a vector to hold the indices of transitions to wake
        wakeTransitions = find(diff(postSleepPeriod) == -postSleepPeriod(1:end-1) & postSleepPeriod(2:end) == 0);
        
        % Calculate the earliest transition to wake after sleep onset, if any
        if ~isempty(wakeTransitions)
            earliestWakeTransition = wakeTransitions(1) * epoch_dur;
        else
            earliestWakeTransition = NaN;
        end

        % Calculate the average sleep depth
        postOnsetSleepDepth = mean(postSleepPeriod); 
        
        % Calculate the fluctuation index for the postSleepPeriod
        fluctuationIndexPost = totalTransitions / totalPostDuration;
        
        % Output post-sleep stats for verification
        fprintf('Post-Sleep Wake Duration (minutes): %g\n', totalWakeDuration);
        fprintf('Post-Sleep NREM1 Duration (minutes): %g\n', totalNREM1DurationPost);
        fprintf('Relative Wake Duration: %g\n', relativeWakeDuration);
        fprintf('Relative NREM1 Duration: %g\n', relativeNREM1Duration);
        fprintf('Earliest Wake Transition (minutes): %g\n', earliestWakeTransition);
        fprintf('Average Sleep Depth: %g\n', postOnsetSleepDepth);
        fprintf('Fluctuation Index (Post-Sleep): %g\n', fluctuationIndexPost);
        
        %% Calculate the variability/peakiness within the S-variable 
        % Fill the missing values in S_variable 
        S_variable = fillmissing(S_variable, "previous");
        S_variable = fillmissing(S_variable, 'next');

        % Calculate STD
        stdValue = std(S_variable);

        % Calculate smoothness (average second derivative of a signal)
        first_derivative = diff(S_variable);
        second_derivative = diff(first_derivative);
        smoothness = mean(abs(second_derivative));

        % Calculate kurtosis 
        kurtosisValue = kurtosis(S_variable);

        % Calculate RMSSD
        differences = diff(S_variable); % calculate successive differences
        squredDifferences = differences.^2; 
        RMSSD = sqrt(mean(squredDifferences));
        
        fprintf('STD: %f\n', stdValue);
        fprintf('Smoothness metric: %f\n', smoothness);
        fprintf('Kurtosis: %f\n', kurtosisValue);
        fprintf('RMSSD: %f\n', RMSSD);


        % Find the row in the table for the current participant and night
        rowIndex = find(strcmp(fitting_result_table.participantName, participantName) & fitting_result_table.nightNumber == n);
        
        % Update the table with both pre-sleep and post-sleep stats
        if ~isempty(rowIndex)
            fitting_result_table(rowIndex, {'EarliestNREM2PreSleep', 'LatestNREM2PreSleep', 'NumberNREM2PreSleep', 'TotalNREM1DurationPreSleep', 'TotalDurationPreSleep', 'ProportionNREM1PreSleep', 'AverageSleepDepthPreSleep','InconsistencyIndex', 'FluctuationIndex', 'RelativeWakeDurationPostSleep', 'RelativeNREM1DurationPostSleep', 'EarliestWakeTransitionPostSleep', 'AverageSleepDepthPostSleep', 'FluctuationIndexPostSleep', 'STD_S', 'Smoothness_S', 'Kurtosis_S', 'RMSSD_S'}) = ...
                table(earliestNREM2BeforeSleep, latestNREM2BeforeSleep, numberPreSleepNREM2, totalNREM1Duration, totalDuration, proportionNREM1, preonsetSleepDepth, inconsistencyIndex, fluctuationIndex, relativeWakeDuration, relativeNREM1Duration, earliestWakeTransition, postOnsetSleepDepth, fluctuationIndexPost, stdValue, smoothness, kurtosisValue, RMSSD);
        end
        
        % Save the updated table as a new file
        newTableName = 'New_fitting_results_sleep_stats_30min.mat';
        saveDir = '/mnt/LongTermStorage/anastasia/Airforce_analysed/Individual Bifurcation Fitting'; % Update this path as needed
        save(fullfile(saveDir, newTableName), 'fitting_result_table');
    end 
end 

% Save the updated table as a new file
newTableName = 'New_fitting_results_sleep_stats_30min.mat';
saveDir = '/mnt/LongTermStorage/anastasia/Airforce_analysed/Individual Bifurcation Fitting'; % Update this path as needed
save(fullfile(saveDir, newTableName), 'fitting_result_table');
