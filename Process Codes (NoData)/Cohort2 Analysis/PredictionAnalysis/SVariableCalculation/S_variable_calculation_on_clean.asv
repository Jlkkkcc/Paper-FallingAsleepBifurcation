%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Author: Anastasia Ilina

% Calculates S-variables from feature-based representations of falling
% asleep EEG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


% Assume participantSubset is your current table with EEG data
Fs = 256; % Sampling rate
python_directory = ''; % Your Python directory
iflinux = true; % Your operating system flag

% Initialize tables to store features for each epoch of each channel
numFeatures = 47;

windowSize = 10; % smooth over 30 seconds

%SelectedChannels = {'F4_A1', 'C4_A1', 'P4_A1', 'O2_A1', 'LOC_A2', 'ROC_A1', 'F3_A2', 'C3_A2', 'P3_A2', 'O1_A2'};
SelectedChannels = {'F4_A1', 'C4_A1', 'P4_A1', 'O2_A1', 'F3_A2', 'C3_A2', 'P3_A2', 'O1_A2'};

sleepOnsetDuration = 60; %in seconds 

saveDirImg = ['/home/anastasia/Dropbox/ISN work/Work on new dataset/Critical Point Prediction/S variable'];
saveDirDerivative = '/home/anastasia/Dropbox/ISN work/Work on new dataset/Critical Point Prediction/Derivative of S Variable';

%saveDirImg = '/home/anastasia/Dropbox/ISN work/Work on new dataset/Full Sleep Onset/S_variable';
%saveDirDerivative = '/home/anastasia/Dropbox/ISN work/Work on new dataset/Full Sleep Onset/Derivative_S_variable';

% Loop through the participant folders to extract the data from the hypnograms 
for i = i:length(participantNames)

    
    participantName = participantNames{i};
    saveDirImgPartInd = [saveDirImg, '/', participantName, '/Local Sleep Onset'];
    saveDirImgPartGlobal = [saveDirImg, '/', participantName, '/Global Sleep Onset'];
    saveDirImgPartIndClean = [saveDirImg, '/', participantName, '/Artifact-Cleaned Local Sleep Onset'];
    saveDirImgPartGlobalClean = [saveDirImg, '/', participantName, '/Artifact-Cleaned Global Sleep Onset'];
    
    saveDirDerivativeInd = [saveDirDerivative, '/', participantName, '/Local Sleep Onset'];
    saveDirDerivativeGlobal = [saveDirDerivative, '/', participantName, '/Global Sleep Onset'];
    saveDirDerivativeIndClean = [saveDirDerivative, '/', participantName, '/Artifact-Cleaned Local Sleep Onset'];
    saveDirDerivativeGlobalClean = [saveDirDerivative, '/', participantName, '/Artifact-Cleaned Global Sleep Onset'];
    
    if ~exist(saveDirImgPartInd, 'dir')
        mkdir(saveDirImgPartInd);
    end

    if ~exist(saveDirImgPartGlobal, 'dir')
        mkdir(saveDirImgPartGlobal);
    end

    if ~exist(saveDirDerivativeInd, 'dir')
        mkdir(saveDirDerivativeInd)
    end 

    if ~exist(saveDirDerivativeGlobal, 'dir')
        mkdir(saveDirDerivativeGlobal)
    end 
    
    if ~exist(saveDirImgPartIndClean, 'dir')
        mkdir(saveDirImgPartIndClean);
    end

    if ~exist(saveDirImgPartGlobalClean, 'dir')
        mkdir(saveDirImgPartGlobalClean);
    end

    if ~exist(saveDirDerivativeIndClean, 'dir')
        mkdir(saveDirDerivativeIndClean)
    end 

    if ~exist(saveDirDerivativeGlobalClean, 'dir')
        mkdir(saveDirDerivativeGlobalClean)
    end

    disp(['S-variable calculation started for Subject No. ',participantName])

    % Construct the path to the .mat file
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
    
    if ~ismember('S_variable', participantSubset.Properties.VariableNames)
        participantSubset.S_variable = cell(height(participantSubset), 1);
    end

    if ~ismember('Clean_S_variable', participantSubset.Properties.VariableNames)
        participantSubset.S_variable = cell(height(participantSubset), 1);
    end

    if ~ismember('Global_S_variable', participantSubset.Properties.VariableNames)
        participantSubset.Global_S_variable = cell(height(participantSubset), 1);
    end

    if ~ismember('Clean_Global_S_variable', participantSubset.Properties.VariableNames)
        participantSubset.Global_S_variable = cell(height(participantSubset), 1);
    end

    if ~ismember('Smoothed_S_variable', participantSubset.Properties.VariableNames)
        participantSubset.Smoothed_S_variable = cell(height(participantSubset), 1);
    end
    
    if ~ismember('Clean_Smoothed_S_variable', participantSubset.Properties.VariableNames)
        participantSubset.Smoothed_S_variable = cell(height(participantSubset), 1);
    end

    if ~ismember('Smoothed_Global_S_variable', participantSubset.Properties.VariableNames)
        participantSubset.Smoothed_Global_S_variable = cell(height(participantSubset), 1);
    end

    if ~ismember('Clean_Smoothed_Global_S_variable', participantSubset.Properties.VariableNames)
        participantSubset.Smoothed_Global_S_variable = cell(height(participantSubset), 1);
    end

    if ~ismember('Derivative_S_variable', participantSubset.Properties.VariableNames)
        participantSubset.Derivative_S_variable = cell(height(participantSubset), 1);
    end

    if ~ismember('Clean_Derivative_S_variable', participantSubset.Properties.VariableNames)
        participantSubset.Derivative_S_variable = cell(height(participantSubset), 1);
    end

    if ~ismember('Global_Derivative_S_variable', participantSubset.Properties.VariableNames)
        participantSubset.Global_Derivative_S_variable = cell(height(participantSubset), 1);
    end

    if ~ismember('Clean_Global_Derivative_S_variable', participantSubset.Properties.VariableNames)
        participantSubset.Global_Derivative_S_variable = cell(height(participantSubset), 1);
    end

    if ~ismember('sleepOnsetIdx', participantSubset.Properties.VariableNames)
        participantSubset.sleepOnsetIdx = NaN(height(participantSubset), 1);
    end

    if ~ismember('timeToSleepOnsetMinutes', participantSubset.Properties.VariableNames)
        participantSubset.timeToSleepOnsetMinutes = cell(height(participantSubset), 1);
    end
 
    cleanFtsAcrossNights = [];
    cleanAvgFtsStorage = {};

    FtsAcrossNights = [];
    avgFtsStorage = {};

     
    %% Calculate mean and std of feature values across all nights for global normalisation later on


    meanFtValuesNight = NaN(1, numFeatures);
    stdFtValuesNight = NaN(1, numFeatures);

    cleanMeanFtValuesNight = NaN(1, numFeatures);
    cleanStdFtValuesNight = NaN(1, numFeatures);

    disp(['Calculating mean and std for participant ...'])

    for n = 1:8 %height(participantSubset)
        if isempty(participantSubset.FallingAsleepFeatures{n})
            continue
        end

        

        FallingAsleepFtsTable = participantSubset.FallingAsleepFeatures{n};
        CleanFallingAsleepFtsTable = participantSubset.CleanFallingAsleepFeatures{n};
        
        % Initialize a matrix to store the features for averaging
        featureMatrixFallingAsleep = NaN(height(participantSubset.FallingAsleepFeatures{n}), numFeatures, length(SelectedChannels));
        cleanFeatureMatrixFallingAsleep = NaN(height(participantSubset.CleanFallingAsleepFeatures{n}), numFeatures, length(SelectedChannels));


        % Loop through the specified EEG channels and accumulate the features
        for ch = 1:length(SelectedChannels)

            channelDataFallingAsleep = participantSubset.FallingAsleepFeatures{n}.(SelectedChannels{ch});
            cleanChannelDataFallingAsleep = participantSubset.CleanFallingAsleepFeatures{n}.(SelectedChannels{ch});
            
            for f = 1:numFeatures

                % Extract the feature values for the current channel
                featureMatrixFallingAsleep(:, f, ch) = cellfun(@(x) x(f), channelDataFallingAsleep);
                cleanFeatureMatrixFallingAsleep(:, f, ch) = cellfun(@(x) x(f), cleanChannelDataFallingAsleep);
    
            end
        end
        
        % Calculate the average features
        averageFeaturesFallingAsleep = mean(featureMatrixFallingAsleep, 3, 'omitmissing');
        cleanAverageFeaturesFallingAsleep = mean(cleanFeatureMatrixFallingAsleep, 3, 'omitmissing');
        avgFtsStorage{n} = averageFeaturesFallingAsleep;
        cleanAvgFtsStorage{n} = cleanAverageFeaturesFallingAsleep;

        % Concatenate the average features for each night
        FtsAcrossNights = [FtsAcrossNights; averageFeaturesFallingAsleep];
        cleanFtsAcrossNights = [cleanFtsAcrossNights; cleanAverageFeaturesFallingAsleep];
    end

    meanFtValuesNight = mean(FtsAcrossNights, 'omitmissing');
    stdFtValuesNight = std(FtsAcrossNights, 'omitmissing');

    cleanMeanFtValuesNight = mean(cleanFtsAcrossNights, 'omitmissing');
    cleanStdFtValuesNight = std(cleanFtsAcrossNights, 'omitmissing');


    %% CALCULATE THE EUCLIDEAN DISTANCE TO SLEEP ONSET (S-VARIABLE) NIGHT-WISE
    
    sleepOnsetAcrossNights = [];
    normalisedFtTrajs = {};

    cleanSleepOnsetAcrossNights = [];
    cleanNormalisedFtTrajs = {};

    max_length = 0;
    
    for n = 1: 8 %height(participantSubset)
        if isempty(participantSubset.FallingAsleepFeatures{n})
            continue
        end
        disp(['Local S-variable calculation started for night ',num2str(n)])
        

        timepoints = participantSubset.FallingAsleepFeatures{n}.Timestamp;
        if length(timepoints) > max_length 
            max_length = length(timepoints);
            max_index = n;
        end 
        

        sleepOnsetSec = (participantSubset.SleepOnsetIndex(n) - 1)*30;
        sleepOnsetIdx = find(timepoints == sleepOnsetSec);
        participantSubset.sleepOnsetIdx(n) = sleepOnsetIdx;
        

        FeatureTrajectoryFallingAsleep = avgFtsStorage{n};
        cleanFeatureTrajectoryFallingAsleep = cleanAvgFtsStorage{n};

        % ADD ARTEFACT CHECK LATER ON
        
        numEpochs = length(FeatureTrajectoryFallingAsleep);
        NormalisedFeatureTrajectoryFallingAsleep = NaN(numEpochs, numFeatures);
        cleanNormalisedFeatureTrajectoryFallingAsleep = NaN(numEpochs, numFeatures);
        distPerNight = NaN(numEpochs, 1);
        cleanDistPerNight = NaN(numEpochs, 1);
        
        % Perform feature normalisation 
        for epch_idx = 1:numEpochs

            NormalisedFeatureTrajectoryFallingAsleep(epch_idx,:) =  (FeatureTrajectoryFallingAsleep(epch_idx, :) - meanFtValuesNight) ./ stdFtValuesNight;
            cleanNormalisedFeatureTrajectoryFallingAsleep(epch_idx,:) =  (cleanFeatureTrajectoryFallingAsleep(epch_idx, :) - cleanMeanFtValuesNight) ./ cleanStdFtValuesNight;
           
        end
        normalisedFtTrajs{n} = NormalisedFeatureTrajectoryFallingAsleep;
        cleanNormalisedFtTrajs{n} = cleanNormalisedFeatureTrajectoryFallingAsleep;

        FtValuesAtOnset = NormalisedFeatureTrajectoryFallingAsleep(sleepOnsetIdx:end,:);
        cleanFtValuesAtOnset = cleanNormalisedFeatureTrajectoryFallingAsleep(sleepOnsetIdx:end,:);
        sleepOnsetAcrossNights = [sleepOnsetAcrossNights; FtValuesAtOnset];
        cleanSleepOnsetAcrossNights = [cleanSleepOnsetAcrossNights; cleanFtValuesAtOnset];
        medianFtsAtOnset = median(FtValuesAtOnset, 1, 'omitmissing');
        cleanMedianFtsAtOnset = median(cleanFtValuesAtOnset, 1, 'omitmissing');

        if sum(isnan(medianFtsAtOnset)) > 0
            continue
        end

        
        
        % Calculate the Euclidean distance
        for epch_idx = 1:numEpochs

            if sum(isnan(NormalisedFeatureTrajectoryFallingAsleep(epch_idx, :))) > 0
                continue
            end 

            distPerNight(epch_idx, :) = pdist2(NormalisedFeatureTrajectoryFallingAsleep(epch_idx, :), medianFtsAtOnset);

        end 

        % Calculate the Euclidean distance
        for epch_idx = 1:numEpochs

            if sum(isnan(cleanNormalisedFeatureTrajectoryFallingAsleep(epch_idx, :))) > 0
                continue
            end 

            cleanDistPerNight(epch_idx, :) = pdist2(cleanNormalisedFeatureTrajectoryFallingAsleep(epch_idx, :), cleanMedianFtsAtOnset);

        end

        
        participantSubset.S_variable{n} = distPerNight;
        smoothedDistPerNight = movmean(distPerNight, windowSize);
        participantSubset.Smoothed_S_variable{n} = smoothedDistPerNight;

        participantSubset.Clean_S_variable{n} = cleanDistPerNight;
        cleanSmoothedDistPerNight = movmean(cleanDistPerNight, windowSize);
        participantSubset.Clean_Smoothed_S_variable{n} = cleanSmoothedDistPerNight;

        % CALCULATE THE DERIVATIVE OF EUCLIDEAN DISTANCE 
        % (if Euclidean distance is not NaN)
        if sum(isnan(distPerNight)) < length(distPerNight)
            derivativeDistPerNight = gradient(smoothedDistPerNight);
            participantSubset.Derivative_S_variable{n} = derivativeDistPerNight;
        end 


        if sum(isnan(cleanDistPerNight)) < length(cleanDistPerNight)
            cleanDerivativeDistPerNight = gradient(cleanSmoothedDistPerNight);
            participantSubset.Clean_Derivative_S_variable{n} = cleanDerivativeDistPerNight;
        end 


        %% CALCULATE TIME TO SLEEP ONSET (IN MINUTES)

        timeToSleepOnsetMinutes = (timepoints - sleepOnsetSec) / 60;
        participantSubset.timeToSleepOnsetMinutes{n} = timeToSleepOnsetMinutes;
        %disp(['Length of timeToSleepOnset is ', num2str(length(timeToSleepOnsetMinutes))])  
        %disp(['Length of distPerNight is ', num2str(length(distPerNight))])

        %% PLOT THE LOCAL S-VARIABLE
        if sum(isnan(distPerNight)) < length(distPerNight)

            

            h = figure; % Open a new figure and store the handle in h
            plot(timeToSleepOnsetMinutes, smoothedDistPerNight)
            
            %set(gca,'FontSize', 20)
            set(gca,'TickDir','out')
            %set(gca,'ticklength',3*get(gca,'ticklength'))
            set(gca,'lineWidth',2)
            xlabel('Time from Sleep Onset (minutes)');
            ylabel('Euclidean Distance from Sleep Onset');
            
            % Construct the title with line breaks
            graphTitle = {'Euclidean distance to sleep onset feature values on this night', ...   
                          ['Participant ', participantName, ' on Night ', num2str(n)], ...
                          ['Took ', num2str(participantSubset.Time2Sleep(n)), ' mins to fall asleep']};
            
            % Set the title
            title(graphTitle, 'Interpreter', 'none');

            ax = gca;
            ylim([0,15])
            line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
            %xlim([0,numEpochs+1])
            
            %xticks([62:100:numEpochs+1])
            %xticklabels([-70:5:10])

            % Construct filenames for saving
            jpegFilename = fullfile(saveDirImgPartInd, ['Participant_', participantName, '_Night', num2str(n), '_Plot.jpg']);
            figFilename = fullfile(saveDirImgPartInd, ['Participant_', participantName, '_Night', num2str(n), '_Plot.fig']);
            
            % Save the plot as JPEG
            saveas(h, jpegFilename);
    
            % Save the figure as .fig file
            savefig(h, figFilename);
    
            close(h); % Close the figure
            
            %% PLOT THE CLEAN S-VARIABLE

            h = figure; % Open a new figure and store the handle in h
            plot(timeToSleepOnsetMinutes, cleanSmoothedDistPerNight)
            
            %set(gca,'FontSize', 20)
            set(gca,'TickDir','out')
            %set(gca,'ticklength',3*get(gca,'ticklength'))
            set(gca,'lineWidth',2)
            xlabel('Time from Sleep Onset (minutes)');
            ylabel('Artifact-Cleaned Euclidean Distance from Sleep Onset');
            
            % Construct the title with line breaks
            graphTitle = {'Euclidean distance to sleep onset feature values on this night', ...   
                          ['Participant ', participantName, ' on Night ', num2str(n)], ...
                          ['Took ', num2str(participantSubset.Time2Sleep(n)), ' mins to fall asleep']};
            
            % Set the title
            title(graphTitle, 'Interpreter', 'none');

            ax = gca;
            ylim([0,15])
            line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
            %xlim([0,numEpochs+1])
            
            %xticks([62:100:numEpochs+1])
            %xticklabels([-70:5:10])

            % Construct filenames for saving
            jpegFilename = fullfile(saveDirImgPartIndClean, ['Participant_', participantName, '_Night', num2str(n), '_Plot.jpg']);
            figFilename = fullfile(saveDirImgPartIndClean, ['Participant_', participantName, '_Night', num2str(n), '_Plot.fig']);
            
            % Save the plot as JPEG
            saveas(h, jpegFilename);
    
            % Save the figure as .fig file
            savefig(h, figFilename);
    
            close(h); % Close the figure

            %% PLOT THE DERIVATIVE OF LOCAL S-VARIABLE

            h = figure; % Open a new figure and store the handle in h
            plot(timeToSleepOnsetMinutes,  movmean(derivativeDistPerNight, 5))

            %set(gca,'FontSize', 20)
            set(gca,'TickDir','out')
            %set(gca,'ticklength',3*get(gca,'ticklength'))
            set(gca,'lineWidth',2)
            xlabel('Time from Sleep Onset (minutes)');
            ylabel('Derivative of Euclidean Distance from Sleep Onset');

             % Construct the title with line breaks
            graphTitle = {'Derivative of Eucl. dist. to sleep onset feature values on this night', ...    
                          ['Participant ', participantName, ' on Night ', num2str(n)], ...
                          ['Took ', num2str(participantSubset.Time2Sleep(n)), ' mins to fall asleep']};
            
            % Set the title
            title(graphTitle, 'Interpreter', 'none');

            ax = gca;
            ylim([-1.5,1.5])
            line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
            xlim([-20,10])
            line(ax.XLim, [0, 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k');
            
            %xticks([62:100:numEpochs+1])
            %xticklabels([-70:5:10])

            % Construct filenames for saving
            jpegFilename = fullfile(saveDirDerivativeInd, ['Participant_', participantName, '_Night', num2str(n), '_Plot.jpg']);
            figFilename = fullfile(saveDirDerivativeInd, ['Participant_', participantName, '_Night', num2str(n), '_Plot.fig']);
            
            % Save the plot as JPEG
            saveas(h, jpegFilename);
    
            % Save the figure as .fig file
            savefig(h, figFilename);
    
            close(h); % Close the figure

            %% PLOT THE ARTIFACT-CLEANED DERIVATIVE OF LOCAL S-VARIABLE

            h = figure; % Open a new figure and store the handle in h
            plot(timeToSleepOnsetMinutes,  movmean(cleanDerivativeDistPerNight, 5))

            %set(gca,'FontSize', 20)
            set(gca,'TickDir','out')
            %set(gca,'ticklength',3*get(gca,'ticklength'))
            set(gca,'lineWidth',2)
            xlabel('Time from Sleep Onset (minutes)');
            ylabel('Derivative of Euclidean Distance from Sleep Onset');

             % Construct the title with line breaks
            graphTitle = {'Art.-Cleaned Derivative of Eucl. dist. to sleep onset feature values on this night', ...    
                          ['Participant ', participantName, ' on Night ', num2str(n)], ...
                          ['Took ', num2str(participantSubset.Time2Sleep(n)), ' mins to fall asleep']};
            
            % Set the title
            title(graphTitle, 'Interpreter', 'none');

            ax = gca;
            ylim([-1.5,1.5])
            line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
            xlim([-20,10])
            line(ax.XLim, [0, 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k');
            
            %xticks([62:100:numEpochs+1])
            %xticklabels([-70:5:10])

            % Construct filenames for saving
            jpegFilename = fullfile(saveDirDerivativeIndClean, ['Participant_', participantName, '_Night', num2str(n), '_Plot.jpg']);
            figFilename = fullfile(saveDirDerivativeIndClean, ['Participant_', participantName, '_Night', num2str(n), '_Plot.fig']);
            
            % Save the plot as JPEG
            saveas(h, jpegFilename);
    
            % Save the figure as .fig file
            savefig(h, figFilename);
    
            close(h); % Close the figure


            
        end      
    
    end
    
    %% CALCULATE EUCLIDEAN DISTANCE TO GLOBAL MEDIAN OF SLEEP ONSET
    

    globalMedianSleepOnset = median(sleepOnsetAcrossNights, 1, 'omitmissing');
    cleanGlobalMedianSleepOnset = median(cleanSleepOnsetAcrossNights, 1, 'omitmissing');

    for n = 1:height(participantSubset)

        if isempty(participantSubset.FallingAsleepFeatures{n})
            continue
        end

        disp(['Global S-variable calculation started for night ',num2str(n)])

    
        if sum(isnan(globalMedianSleepOnset)) > 0
            continue
        end
        
        FtTrajectory = normalisedFtTrajs{n};
        cleanFtTrajectory = cleanNormalisedFtTrajs{n};
        timepoints = participantSubset.FallingAsleepFeatures{n}.Timestamp;
        sleepOnsetSec = (participantSubset.SleepOnsetIndex(n) - 1)*30;

        
        numEpochs = length(FtTrajectory);
        distPerNight = NaN(numEpochs, 1);
        cleanDistPerNight = NaN(numEpochs, 1);
        
        % Calculate the Euclidean distance
        for epch_idx = 1:numEpochs

            if sum(isnan(FtTrajectory(epch_idx, :))) > 0
                continue
            end 

            distPerNight(epch_idx, :) = pdist2(FtTrajectory(epch_idx, :), globalMedianSleepOnset);

        end 

        % Calculate the Artifact-Cleaned Euclidean distance
        for epch_idx = 1:numEpochs

            if sum(isnan(cleanFtTrajectory(epch_idx, :))) > 0
                continue
            end 

            cleanDistPerNight(epch_idx, :) = pdist2(cleanFtTrajectory(epch_idx, :), cleanGlobalMedianSleepOnset);

        end

        
        participantSubset.Global_S_variable{n} = distPerNight;
        smoothedDistPerNight = movmean(distPerNight, windowSize);
        participantSubset.Smoothed_Global_S_variable{n} = smoothedDistPerNight;

        participantSubset.Clean_Global_S_variable{n} = cleanDistPerNight;
        cleanSmoothedDistPerNight = movmean(cleanDistPerNight, windowSize);
        participantSubset.Clean_Smoothed_Global_S_variable{n} = cleanSmoothedDistPerNight;
        
        %% CALCULATE THE DERIVATIVE OF GLOBAL S-VARIABLE
        % (if Euclidean distance is not NaN)
        if sum(isnan(distPerNight)) < length(distPerNight)
            derivativeDistPerNight = gradient(smoothedDistPerNight);
            participantSubset.Global_Derivative_S_variable{n} = derivativeDistPerNight;
        end 

        %% CALCULATE THE ARTIFACT-CLEANED DERIVATIVE OF GLOBAL S-VARIABLE
        % (if Euclidean distance is not NaN)
        if sum(isnan(cleanDistPerNight)) < length(cleanDistPerNight)
            cleanDerivativeDistPerNight = gradient(cleanSmoothedDistPerNight);
            participantSubset.Clean_Global_Derivative_S_variable{n} = cleanDerivativeDistPerNight;
        end 

        %% CALCULATE TIME TO SLEEP ONSET
        timeToSleepOnsetMinutes = (timepoints - sleepOnsetSec) / 60;
        %disp(['Length of timeToSleepOnset is ', num2str(length(timeToSleepOnsetMinutes))])  
        %disp(['Length of distPerNight is ', num2str(length(distPerNight))])

        %% PLOT GLOBAL S-VARIABLE FOR THE NIGHT
        if sum(isnan(distPerNight)) < length(distPerNight)
            h = figure; % Open a new figure and store the handle in h
            plot(timeToSleepOnsetMinutes,  smoothedDistPerNight)
            
            %set(gca,'FontSize', 20)
            set(gca,'TickDir','out')
            %set(gca,'ticklength',3*get(gca,'ticklength'))
            set(gca,'lineWidth',2)
            xlabel('Time from Sleep Onset (minutes)');
            ylabel('Euclidean Distance from Global Sleep Onset');
            
            % Construct the title with line breaks
            graphTitle = {'Euclidean distance to global sleep onset feature values', ...    
                          ['Participant ', participantName, ' on Night ', num2str(n)], ...
                          ['Took ', num2str(participantSubset.Time2Sleep(n)), ' mins to fall asleep']};
            
            % Set the title
            title(graphTitle, 'Interpreter', 'none');

            ax = gca;
            ylim([0,15])
            line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
            %xlim([0,numEpochs+1])
            
            %xticks([62:100:numEpochs+1])
            %xticklabels([-70:5:10])

            % Construct filenames for saving
            jpegFilename = fullfile(saveDirImgPartGlobal, ['Participant_', participantName, '_Night', num2str(n), '_Plot.jpg']);
            figFilename = fullfile(saveDirImgPartGlobal, ['Participant_', participantName, '_Night', num2str(n), '_Plot.fig']);
            
            % Save the plot as JPEG
            saveas(h, jpegFilename);
    
            % Save the figure as .fig file
            savefig(h, figFilename);
    
            close(h); % Close the figure

            %% PLOT THE ARTIFACT-CLEANED GLOBAL S-VARIABLE FOR THE NIGHT
            h = figure; % Open a new figure and store the handle in h
            plot(timeToSleepOnsetMinutes,  cleanSmoothedDistPerNight)
            
            %set(gca,'FontSize', 20)
            set(gca,'TickDir','out')
            %set(gca,'ticklength',3*get(gca,'ticklength'))
            set(gca,'lineWidth',2)
            xlabel('Time from Sleep Onset (minutes)');
            ylabel('Euclidean Distance from Global Sleep Onset');
            
            % Construct the title with line breaks
            graphTitle = {'Artifact-Cleaned Euclidean distance to global sleep onset feature values', ...    
                          ['Participant ', participantName, ' on Night ', num2str(n)], ...
                          ['Took ', num2str(participantSubset.Time2Sleep(n)), ' mins to fall asleep']};
            
            % Set the title
            title(graphTitle, 'Interpreter', 'none');

            ax = gca;
            ylim([0,15])
            line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
            %xlim([0,numEpochs+1])
            
            %xticks([62:100:numEpochs+1])
            %xticklabels([-70:5:10])

            % Construct filenames for saving
            jpegFilename = fullfile(saveDirImgPartGlobalClean, ['Participant_', participantName, '_Night', num2str(n), '_Plot.jpg']);
            figFilename = fullfile(saveDirImgPartGlobalClean, ['Participant_', participantName, '_Night', num2str(n), '_Plot.fig']);
            
            % Save the plot as JPEG
            saveas(h, jpegFilename);
    
            % Save the figure as .fig file
            savefig(h, figFilename);
    
            close(h); % Close the figure

            %% PLOT THE DERIVATIVE OF GLOBAL S-VARIABLE FOR THE NIGHT

            h = figure; % Open a new figure and store the handle in h
            plot(timeToSleepOnsetMinutes,  movmean(derivativeDistPerNight, 5))

            %set(gca,'FontSize', 20)
            set(gca,'TickDir','out')
            %set(gca,'ticklength',3*get(gca,'ticklength'))
            set(gca,'lineWidth',2)
            xlabel('Time from Sleep Onset (minutes)');
            ylabel('Derivative of Euclidean Distance from Sleep Onset');

             % Construct the title with line breaks
            graphTitle = {'Derivative of Eucl. dist. to global sleep onset feature values', ...    
                          ['Participant ', participantName, ' on Night ', num2str(n)], ...
                          ['Took ', num2str(participantSubset.Time2Sleep(n)), ' mins to fall asleep']};
            
            % Set the title
            title(graphTitle, 'Interpreter', 'none');

            ax = gca;
            ylim([-1.5,1.5])
            line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
            line(ax.XLim, [0, 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k');
            xlim([-20,10])
            
            %xticks([62:100:numEpochs+1])
            %xticklabels([-70:5:10])

            % Construct filenames for saving
            jpegFilename = fullfile(saveDirDerivativeGlobal, ['Participant_', participantName, '_Night', num2str(n), '_Plot.jpg']);
            figFilename = fullfile(saveDirDerivativeGlobal, ['Participant_', participantName, '_Night', num2str(n), '_Plot.fig']);
            
            % Save the plot as JPEG
            saveas(h, jpegFilename);
    
            % Save the figure as .fig file
            savefig(h, figFilename);
    
            close(h); % Close the figure

            %% PLOT THE ARTIFACT-CLEANED DERIVATIVE OF GLOBAL S-VARIABLE FOR THE NIGHT

            h = figure; % Open a new figure and store the handle in h
            plot(timeToSleepOnsetMinutes,  movmean(cleanDerivativeDistPerNight, 5))

            %set(gca,'FontSize', 20)
            set(gca,'TickDir','out')
            %set(gca,'ticklength',3*get(gca,'ticklength'))
            set(gca,'lineWidth',2)
            xlabel('Time from Sleep Onset (minutes)');
            ylabel('Derivative of Euclidean Distance from Sleep Onset');

             % Construct the title with line breaks
            graphTitle = {'Art.-Cleaned Derivative of Eucl. dist. to global sleep onset feature values', ...    
                          ['Participant ', participantName, ' on Night ', num2str(n)], ...
                          ['Took ', num2str(participantSubset.Time2Sleep(n)), ' mins to fall asleep']};
            
            % Set the title
            title(graphTitle, 'Interpreter', 'none');

            ax = gca;
            ylim([-1.5,1.5])
            line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
            line(ax.XLim, [0, 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k');
            xlim([-20,10])
            
            %xticks([62:100:numEpochs+1])
            %xticklabels([-70:5:10])

            % Construct filenames for saving
            jpegFilename = fullfile(saveDirDerivativeGlobalClean, ['Participant_', participantName, '_Night', num2str(n), '_Plot.jpg']);
            figFilename = fullfile(saveDirDerivativeGlobalClean, ['Participant_', participantName, '_Night', num2str(n), '_Plot.fig']);
            
            % Save the plot as JPEG
            saveas(h, jpegFilename);
    
            % Save the figure as .fig file
            savefig(h, figFilename);
    
            close(h); % Close the figure
        end 
        
        
    end
    % Save the participant-specific subset table using MAT-file version 7.3
    savePathParticipantSubset = fullfile('/mnt/LongTermStorage/anastasia/Airforce_analysed/Participants_Epoched_EEG/', ['ArtifactFtEEGFallingAsleep_', participantName, '.mat']);

    % Use the '-v7.3' flag to save larger data
    save(savePathParticipantSubset, 'participantSubset', '-v7.3');
    
    %% POPULATE MATRICES FOR GLOBAL AND LOCAL S-VARIABLES AND ITS DERIVARIVES
    % FOR MEAN AND STD CALCLATION

    allSVariables = NaN(n, max_length);
    allGlobalSVariables = NaN(n, max_length);
    allDerivativesSLocal = NaN(n, max_length);
    allDerivativesSGlobal = NaN(n, max_length);

    allSVariablesClean = NaN(n, max_length);
    allGlobalSVariablesClean = NaN(n, max_length);
    allDerivativesSLocalClean = NaN(n, max_length);
    allDerivativesSGlobalClean = NaN(n, max_length);

    for n = 1:height(participantSubset)
        S_variable = participantSubset.S_variable{n};
        DerivativeLocal = participantSubset.Derivative_S_variable{n};
        Clean_S_variable = participantSubset.Clean_S_variable{n};
        CleanDerivativeLocal = participantSubset.Clean_Derivative_S_variable{n};
        numEpochsLocal = length(S_variable);

        Global_S_variable = participantSubset.Global_S_variable{n};
        DerivativeGlobal = participantSubset.Global_Derivative_S_variable{n};
        Clean_Global_S_variable = participantSubset.Clean_Global_S_variable{n};
        CleanDerivativeGlobal = participantSubset.Clean_Global_Derivative_S_variable{n};
        numEpochsGlobal = length(Global_S_variable);
        
        allSVariables(n, 1+max_length-numEpochsLocal:end) = S_variable;
        allDerivativesSLocal(n, 1+max_length-numEpochsLocal:end) = DerivativeLocal;
        allGlobalSVariables(n, 1+max_length-numEpochsGlobal:end) = Global_S_variable;
        allDerivativesSGlobal(n, 1+max_length-numEpochsGlobal:end) = DerivativeGlobal;

        allSVariablesClean(n, 1+max_length-numEpochsLocal:end) = Clean_S_variable;
        allDerivativesSLocalClean(n, 1+max_length-numEpochsLocal:end) = CleanDerivativeLocal;
        allGlobalSVariablesClean(n, 1+max_length-numEpochsGlobal:end) = Clean_Global_S_variable;
        allDerivativesSGlobalClean(n, 1+max_length-numEpochsGlobal:end) = CleanDerivativeGlobal;
    end
    
    %% MEAN OF LOCAL S-VARIABLE AND ITS DERIVATIVE ACROSS NIGHTS
    disp('Mean Local S-variable calculation ...')


    meanLocalSVariable = mean(allSVariables,1, 'omitmissing');
    stdLocalSVariabel = std(allSVariables,1, 'omitmissing');
    meanLocalSDerivative = mean(allDerivativesSLocal, 1, 'omitmissing');
    stdLocalSDerivative = std(allDerivativesSLocal, 1, 'omitmissing');

    cleanMeanLocalSVariable = mean(allSVariables,1, 'omitmissing');
    cleanStdLocalSVariabel = std(allSVariables,1, 'omitmissing');
    cleanMeanLocalSDerivative = mean(allDerivativesSLocal, 1, 'omitmissing');
    cleanStdLocalSDerivative = std(allDerivativesSLocal, 1, 'omitmissing');

    timepoints_max = participantSubset.FallingAsleepFeatures{max_index}.Timestamp;
    sleepOnsetSec_max = (participantSubset.SleepOnsetIndex(max_index) - 1)*30;
    MaxTimeToSleepOnsetMinutes = (timepoints_max - sleepOnsetSec_max) / 60;
    
    %% PLOT MEAN AND STD OF LOCAL S-VARIABLE
    h = figure; % Open a new figure and store the handle in h
    shadedErrorBar(MaxTimeToSleepOnsetMinutes,  meanLocalSVariable, stdLocalSVariabel, 'lineProps', {'-k', 'lineWidth', 0.5})
    
    %set(gca,'FontSize', 20)
    set(gca,'TickDir','out')
    %set(gca,'ticklength',3*get(gca,'ticklength'))
    set(gca,'lineWidth',2)
    xlabel('Time from Sleep Onset (minutes)');
    ylabel('Mean Euclidean Distance from Sleep Onset');
    
    % Construct the title with line breaks
    graphTitle = {['Participant ', participantName, ' across all Nights ', num2str(n)]};
                  
    
    % Set the title
    title(graphTitle, 'Interpreter', 'none');

    ax = gca;
    ylim([0,15])
    line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
    xlim([-60, 10])
    
    %xticks([62:100:numEpochs+1])
    %xticklabels([-70:5:10])

    % Construct filenames for saving
    jpegFilename = fullfile(saveDirImgPartInd, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.jpg']);
    figFilename = fullfile(saveDirImgPartInd, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.fig']);
    
    % Save the plot as JPEG
    saveas(h, jpegFilename);

    % Save the figure as .fig file
    savefig(h, figFilename);

    close(h); % Close the figure

     %% PLOT MEAN AND STD OF ARTIFACT-CLEANED LOCAL S-VARIABLE
    h = figure; % Open a new figure and store the handle in h
    shadedErrorBar(MaxTimeToSleepOnsetMinutes,  cleanMeanLocalSVariable, cleanStdLocalSVariabel, 'lineProps', {'-k', 'lineWidth', 0.5})
    
    %set(gca,'FontSize', 20)
    set(gca,'TickDir','out')
    %set(gca,'ticklength',3*get(gca,'ticklength'))
    set(gca,'lineWidth',2)
    xlabel('Time from Sleep Onset (minutes)');
    ylabel('Artifact-Cleaned Mean Euclidean Distance from Sleep Onset');
    
    % Construct the title with line breaks
    graphTitle = {['Participant ', participantName, ' across all Nights ', num2str(n)]};
                  
    
    % Set the title
    title(graphTitle, 'Interpreter', 'none');

    ax = gca;
    ylim([0,15])
    line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
    xlim([-60, 10])
    
    %xticks([62:100:numEpochs+1])
    %xticklabels([-70:5:10])

    % Construct filenames for saving
    jpegFilename = fullfile(saveDirImgPartIndClean, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.jpg']);
    figFilename = fullfile(saveDirImgPartIndClean, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.fig']);
    
    % Save the plot as JPEG
    saveas(h, jpegFilename);

    % Save the figure as .fig file
    savefig(h, figFilename);

    close(h); % Close the figure

    %% PLOT AVERAGE DERIVATIVE OF LOCAL S-VARIABLE ACROSS THE NIGHTS 

    h = figure; % Open a new figure and store the handle in h
    shadedErrorBar(MaxTimeToSleepOnsetMinutes,  meanLocalSDerivative, stdLocalSDerivative, 'lineProps', {'-k', 'lineWidth', 0.5})
    
    %set(gca,'FontSize', 20)
    set(gca,'TickDir','out')
    %set(gca,'ticklength',3*get(gca,'ticklength'))
    set(gca,'lineWidth',2)
    xlabel('Time from Sleep Onset (minutes)');
    ylabel('Mean Derivative of Euclidean Distance from Sleep Onset');
    
    % Construct the title with line breaks
    graphTitle = {'Mean Derivative of Eucl. Dist. from Sleep Onset feature values on that night'
                   ['Participant ', participantName, ' across all Nights ', num2str(n)]};
                  
    
    % Set the title
    title(graphTitle, 'Interpreter', 'none');

    ax = gca;
    ylim([-1.5,1.5])
    line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
    line(ax.XLim, [0, 0], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'b');
    xlim([-20,10])
    
    %xticks([62:100:numEpochs+1])
    %xticklabels([-70:5:10])

    % Construct filenames for saving
    jpegFilename = fullfile(saveDirDerivativeInd, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.jpg']);
    figFilename = fullfile(saveDirDerivativeInd, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.fig']);
    
    % Save the plot as JPEG
    saveas(h, jpegFilename);

    % Save the figure as .fig file
    savefig(h, figFilename);

    close(h); % Close the figure

    %% PLOT AVERAGE ARTIFACT-CLEANED DERIVATIVE OF LOCAL S-VARIABLE ACROSS THE NIGHTS 

    h = figure; % Open a new figure and store the handle in h
    shadedErrorBar(MaxTimeToSleepOnsetMinutes,  cleanMeanLocalSDerivative, cleanStdLocalSDerivative, 'lineProps', {'-k', 'lineWidth', 0.5})
    
    %set(gca,'FontSize', 20)
    set(gca,'TickDir','out')
    %set(gca,'ticklength',3*get(gca,'ticklength'))
    set(gca,'lineWidth',2)
    xlabel('Time from Sleep Onset (minutes)');
    ylabel('Mean Derivative of Euclidean Distance from Sleep Onset');
    
    % Construct the title with line breaks
    graphTitle = {'Mean Art-Cleaned Derivative of Eucl. Dist. from Sleep Onset feature values on that night'
                   ['Participant ', participantName, ' across all Nights ', num2str(n)]};
                  
    
    % Set the title
    title(graphTitle, 'Interpreter', 'none');

    ax = gca;
    ylim([-1.5,1.5])
    line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
    line(ax.XLim, [0, 0], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'b');
    xlim([-20,10])
    
    %xticks([62:100:numEpochs+1])
    %xticklabels([-70:5:10])

    % Construct filenames for saving
    jpegFilename = fullfile(saveDirDerivativeIndClean, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.jpg']);
    figFilename = fullfile(saveDirDerivativeIndClean, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.fig']);
    
    % Save the plot as JPEG
    saveas(h, jpegFilename);

    % Save the figure as .fig file
    savefig(h, figFilename);

    close(h); % Close the figure

    %% MEAN OF GLOBAL S-VARIABLE AND ITS DERIVATIVE
     disp('Mean Global S-variable calculation ...')
    meanGlobalSVariable = mean(allGlobalSVariables,1, 'omitmissing');
    stdGlobalSVariabel = std(allGlobalSVariables,1, 'omitmissing');
    meanGlobalSDerivative = mean(allDerivativesSGlobal, 1, 'omitmissing');
    stdGlobalSDerivative = std(allDerivativesSGlobal, 1, 'omitmissing');
    
    cleanMeanGlobalSVariable = mean(allGlobalSVariablesClean,1, 'omitmissing');
    cleanStdGlobalSVariabel = std(allGlobalSVariablesClean,1, 'omitmissing');
    cleanMeanGlobalSDerivative = mean(allDerivativesSGlobalClean, 1, 'omitmissing');
    cleanStdGlobalSDerivative = std(allDerivativesSGlobalClean, 1, 'omitmissing');

    h = figure; % Open a new figure and store the handle in h
    shadedErrorBar(MaxTimeToSleepOnsetMinutes,  meanGlobalSVariable, stdGlobalSVariabel, 'lineProps', {'-k', 'lineWidth', 0.5})
    
    %set(gca,'FontSize', 20)
    set(gca,'TickDir','out')
    %set(gca,'ticklength',3*get(gca,'ticklength'))
    set(gca,'lineWidth',2)
    xlabel('Time from Sleep Onset (minutes)');
    ylabel('Mean Euclidean Distance from Global Sleep Onset');
    
    % Construct the title with line breaks
    graphTitle = {['Participant ', participantName, ' across all Nights ', num2str(n)]};
                  
    
    % Set the title
    title(graphTitle, 'Interpreter', 'none');

    ax = gca;
    ylim([0,15])
    line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
    xlim([-60, 10])
    
    %xticks([62:100:numEpochs+1])
    %xticklabels([-70:5:10])

    % Construct filenames for saving
    jpegFilename = fullfile(saveDirImgPartGlobal, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.jpg']);
    figFilename = fullfile(saveDirImgPartGlobal, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.fig']);
    
    % Save the plot as JPEG
    saveas(h, jpegFilename);

    % Save the figure as .fig file
    savefig(h, figFilename);

    close(h); % Close the figure
    
    % Artefact-cleaned global S-variable
    h = figure; % Open a new figure and store the handle in h
    shadedErrorBar(MaxTimeToSleepOnsetMinutes,  cleanMeanGlobalSVariable, cleanStdGlobalSVariabel, 'lineProps', {'-k', 'lineWidth', 0.5})
    
    %set(gca,'FontSize', 20)
    set(gca,'TickDir','out')
    %set(gca,'ticklength',3*get(gca,'ticklength'))
    set(gca,'lineWidth',2)
    xlabel('Time from Sleep Onset (minutes)');
    ylabel('Artifact-Cleaned Mean Euclidean Distance from Global Sleep Onset');
    
    % Construct the title with line breaks
    graphTitle = {['Participant ', participantName, ' across all Nights ', num2str(n)]};
                  
    
    % Set the title
    title(graphTitle, 'Interpreter', 'none');

    ax = gca;
    ylim([0,15])
    line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
    xlim([-60, 10])
    
    %xticks([62:100:numEpochs+1])
    %xticklabels([-70:5:10])

    % Construct filenames for saving
    jpegFilename = fullfile(saveDirImgPartGlobalClean, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.jpg']);
    figFilename = fullfile(saveDirImgPartGlobalClean, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.fig']);
    
    % Save the plot as JPEG
    saveas(h, jpegFilename);

    % Save the figure as .fig file
    savefig(h, figFilename);

    close(h); % Close the figure

     % Average derivative across nights 

    h = figure; % Open a new figure and store the handle in h
    shadedErrorBar(MaxTimeToSleepOnsetMinutes,  meanGlobalSDerivative, stdGlobalSDerivative, 'lineProps', {'-k', 'lineWidth', 0.5})
    
    %set(gca,'FontSize', 20)
    set(gca,'TickDir','out')
    %set(gca,'ticklength',3*get(gca,'ticklength'))
    set(gca,'lineWidth',2)
    xlabel('Time from Sleep Onset (minutes)');
    ylabel('Mean Derivative of Euclidean Distance from Sleep Onset');
    
    % Construct the title with line breaks
    graphTitle = {'Mean Derivative of Eucl. Dist. from Global Sleep Onset feature values'
                   ['Participant ', participantName, ' across all Nights ', num2str(n)]};
                  
    
    % Set the title
    title(graphTitle, 'Interpreter', 'none');

    ax = gca;
    ylim([-1.5,1.5])
    line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
    line(ax.XLim, [0, 0], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'b');
    xlim([-20,10])
    
    %xticks([62:100:numEpochs+1])
    %xticklabels([-70:5:10])

    % Construct filenames for saving
    jpegFilename = fullfile(saveDirDerivativeGlobal, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.jpg']);
    figFilename = fullfile(saveDirDerivativeGlobal, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.fig']);
    
    % Save the plot as JPEG
    saveas(h, jpegFilename);

    % Save the figure as .fig file
    savefig(h, figFilename);

    close(h); % Close the figure

     % Average derivative across nights 

    h = figure; % Open a new figure and store the handle in h
    shadedErrorBar(MaxTimeToSleepOnsetMinutes,  cleanMeanGlobalSDerivative, cleanStdGlobalSDerivative, 'lineProps', {'-k', 'lineWidth', 0.5})
    
    %set(gca,'FontSize', 20)
    set(gca,'TickDir','out')
    %set(gca,'ticklength',3*get(gca,'ticklength'))
    set(gca,'lineWidth',2)
    xlabel('Time from Sleep Onset (minutes)');
    ylabel('Mean Derivative of Euclidean Distance from Sleep Onset');
    
    % Construct the title with line breaks
    graphTitle = {'Mean Art-Cleaned Derivative of Eucl. Dist. from Global Sleep Onset feature values'
                   ['Participant ', participantName, ' across all Nights ', num2str(n)]};
                  
    
    % Set the title
    title(graphTitle, 'Interpreter', 'none');

    ax = gca;
    ylim([-1.5,1.5])
    line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
    line(ax.XLim, [0, 0], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'b');
    xlim([-20,10])
    
    %xticks([62:100:numEpochs+1])
    %xticklabels([-70:5:10])

    % Construct filenames for saving
    jpegFilename = fullfile(saveDirDerivativeGlobalClean, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.jpg']);
    figFilename = fullfile(saveDirDerivativeGlobalClean, ['Participant_', participantName, '_AllNights', num2str(n), '_Plot.fig']);
    
    % Save the plot as JPEG
    saveas(h, jpegFilename);

    % Save the figure as .fig file
    savefig(h, figFilename);

    close(h); % Close the figure
    
end 


