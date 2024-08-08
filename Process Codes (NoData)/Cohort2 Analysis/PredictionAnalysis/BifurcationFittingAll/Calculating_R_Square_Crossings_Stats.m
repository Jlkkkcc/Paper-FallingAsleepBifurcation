%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Author: Anastasia Ilina

% Fits a bifuraction to true S-variable from all participants and gathers
% statistics about the bifuraction fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fixing directories

load excluded_participants.mat
%load participantSubset.mat
% Add necessary directories
%addpath('/home/anastasia/Dropbox/ISN work/Work on new dataset/Mean S-variable for prediction/Functions')
%addpath('/home/anastasia/Dropbox/ISN work/Work on new dataset')
addpath('/mnt/LongTermStorage/Airforce_Surrey_psa/af0165')
%addpath ('/home/anastasia/Dropbox/ISN work/Whole Code Pipelines/My code/edf-viewer-master')
addpath('/mnt/LongTermStorage/Airforce_Surrey_psa/HYPNOGRAMS')
addpath('/mnt/LongTermStorage/anastasia/Airforce_analysed')
%addpath('/home/anastasia/Dropbox/Model fitting')
addpath('/home/anastasia/Dropbox/ISN work/Work on new dataset/Critical Point Prediction Corrected/')


%% Get a list of all of the participants in the directory
% change directory to the folder with the hypnograms 
cd('/mnt/LongTermStorage/Airforce_Surrey_psa')

% Get a list of all files and folders in the current directory
files = dir;

% Filter out all the files (keep only folders)
folders = files([files.isdir]);

% Remove the '.' and '..' directories
folders = folders(~ismember({folders.name}, {'.', '..', 'HYPNOGRAMS'}));

% Extract the folder names
participantNames = {folders.name};

%% Further initalisation
num_nights = 8;

% Turn off the warnings to declutter the command line
warning('off', 'all')
smoothingWindow = 20;

% Define the constants for the script
epc_len = 6; % epoch length in seconds
overlap = 3; % overlap in seconds
maxTimePreSleep = 30; % maximal pre-sleep time in minutes
maxTimePostSleep = 10; % maximal post-sleep time in minutes
SelectedChannels = {'F4_A1', 'C4_A1', 'P4_A1', 'O2_A1', 'F3_A2', 'C3_A2', 'P3_A2', 'O1_A2'};

smoothingMethod = 'retrospective_median';
smoothingWindow = 20;
smoothingOrder = 4;

numFts = 47;
step_size_sec = epc_len - overlap;

numEpochsPreSleep = (maxTimePreSleep*60)/step_size_sec -1;
numEpochsPostSleep = (maxTimePostSleep * 60)/step_size_sec - 1;


% Define the directories to save the resultant plots and variables
saveDir = '/home/anastasia/Dropbox/ISN work/Work on new dataset/Critical Point Prediction Corrected/';

% Construct path to the Results Table
filePath = fullfile('/home/anastasia/Dropbox/ISN work/Work on new dataset/Critical Point Prediction Corrected/FittingStatsTable.mat');

if exist(filePath, 'file')
    % Load the table from the file
    loadedData = load(filePath);
    fittingStatsTable = loadedData.fittingStatsTable;
else

    fittingStatsTable = table('Size', [0, 19], ...
        'VariableTypes', {'string', 'double', 'double', 'cell', 'cell', 'cell', 'double', 'double', 'cell', 'cell', 'cell', 'double',  'cell', 'logical', 'double', 'double', 'cell', 'double', 'double'}, ...
        'VariableNames', {'Patient_ID', 'Max_Number_Nights_Available', 'Night_Number', 'Local_S_Variable', 'Hypnogram', 'UpdatedTimeVec', 'R_square', 'MSE', 'Optimal_Params', ...
        'Fitted_Bifurcation', 'Control_Parameter', 'SVariableCritical', 'C_3_Solutions', 'ifNoCrossingPoints', 'FinalCrossingTime', 'FalseAlarmRate', 'CrossingDetails', 'LastPositiveCrossing', 'FalseAlarmN1N2Proportion'});

    % Optionally, save the  new table to the file
    save(filePath, 'fittingStatsTable');

end

saveDir = '/home/anastasia/Dropbox/ISN work/Work on new dataset/Critical Point Prediction Corrected/';

saveDirAdjustedBifurcation = [saveDir, 'FittingBifurcations'];
if ~exist(saveDirAdjustedBifurcation, 'dir')
   mkdir(saveDirAdjustedBifurcation);
end

% 
% if restart == 1
%     averageCosineSimilarityScoresAll = [];
%     cosineSimilarityAll = [];
%     num_night_vec_all = [];
%     all_nums = [];
%     average_cosine_scores_mat = NaN(length(participantNames), num_nights);
%     excluded_participants = {};
% else
%     %excluded_participants = {participantNames{8}, participantNames{5}};
%     load(fullfile([saveDir, 'averageCosineSimilarityScoresAll.mat']));
%     load(fullfile([saveDir, 'cosineSimilarityAll.mat']));
%     load(fullfile([saveDir, 'num_night_vec_all.mat']));
%     load(fullfile([saveDir, 'all_nums.mat']));
%     load(fullfile([saveDir, 'average_cosine_scores_mat.mat']));
%     load(fullfile([saveDir, 'excluded_participants.mat']));
% 
% end

%% Run calculations for all participants 

table_idx = 0;
for part = 1: length(participantNames)
    

    participantName = participantNames{part};
    disp(['Analysing participant ', participantName]);


    if part == 8 || part == 5
        excluded_participants = [excluded_participants, participantName];
        continue
    end

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

    % Identifiy which nights are not empty and have acceptable sleep
    % quality after sleep onsest
    nights = 1:8;
    non_empty_mask = ones(num_nights,1);

    for n_n = 1:num_nights

        allNightHypnogram = participantSubset.CleanSleepStagesInt{n_n};
        sleepOnsetIndex = participantSubset.SleepOnsetIndex(n_n);

        if isempty(participantSubset.CleanSleepStagesInt{n_n}) || (participantSubset.SleepOnsetIndex(n_n) + 10*2-1) > length(allNightHypnogram)
            non_empty_mask(n_n) = 0;
            continue
        end


        postSleepPeriod = allNightHypnogram(sleepOnsetIndex:sleepOnsetIndex + 10*2-1);
        preSleepPeriod = allNightHypnogram(1:sleepOnsetIndex -1);
        averagePostSleepDepth = mean(postSleepPeriod);
        numAwakeEpochs = sum(postSleepPeriod == 0);
        REMbeforeSO = sum(preSleepPeriod == 5);


        if averagePostSleepDepth < 1.5 || numAwakeEpochs > 1 || isempty(participantSubset.FallingAsleepFeatures{n_n}) || REMbeforeSO > 0
            non_empty_mask(n_n) = 0;
        end
    end

    non_empty_nights = nights(logical(non_empty_mask(1:8)));

    if sum(non_empty_mask) < 5
        excluded_participants = [excluded_participants, participantName];
        continue
    end


    fittingStatsTableRow = table('Size', [1, 19], ...
        'VariableTypes', {'string', 'double', 'double', 'cell', 'cell', 'cell', 'double', 'double', 'cell', 'cell', 'cell', 'double',  'cell', 'logical', 'double', 'double', 'cell', 'double', 'double'}, ...
        'VariableNames', {'Patient_ID', 'Max_Number_Nights_Available', 'Night_Number', 'Local_S_Variable', 'Hypnogram', 'UpdatedTimeVec', 'R_square', 'MSE', 'Optimal_Params', ...
        'Fitted_Bifurcation', 'Control_Parameter', 'SVariableCritical', 'C_3_Solutions', 'ifNoCrossingPoints', 'FinalCrossingTime', 'FalseAlarmRate', 'CrossingDetails', 'LastPositiveCrossing', 'FalseAlarmN1N2Proportion'});

    for i = 1:length(non_empty_nights)

        table_idx = table_idx+1;

        curIdx = non_empty_nights(i);
         disp(['Analysing Night ', num2str(curIdx)]);

        %% Adjust the bifurcation fittings to have proper critical points

        distPerNight1 = participantSubset.S_variable{curIdx};
        uncut_xx_smoothed = smoothdata(distPerNight1, 'movmedian', [smoothingWindow, 0], 'omitmissing');
        uncut_xx_smoothed  = uncut_xx_smoothed - min(uncut_xx_smoothed);
        uncut_xx_smoothed = uncut_xx_smoothed(1:length(uncut_xx_smoothed) - 9);
        %uncut_tvec = linspace(-30, 10, length(uncut_xx_smoothed));
        uncut_tvec =  (-30*60 + 6):3:10*60;
        uncut_tvec = uncut_tvec / 60;

     
        tvec = uncut_tvec(max(1, length(uncut_tvec) - length(uncut_xx_smoothed)+1): end); 
        xx_smoothed = (uncut_xx_smoothed(max(1,  length(uncut_xx_smoothed) - length(tvec) +1): end))';
        



        saveDirParticipant = [saveDirAdjustedBifurcation, participantName, '/' ];
        if ~exist(saveDirParticipant, 'dir')
            mkdir(saveDirParticipant);
        end

        %% Plot the Orig Bifurcation
        % bif = uniqueTable.Fitted_Bifurcation_Model{curIdx};
        % bif = bif(1:length(bif) - 9);
        % c = uniqueTable.Control_Parameter{curIdx};
        % c = c(1:length(c) - 9);
        %
        % dd_orig = [bif,c];

        % init_params = uniqueTable.Optimal_Params{curIdx};
        % if isnan(uncut_xx_smoothed(1))
        %     plot_bifurcation(dd_orig(firstNonNaNIndex:end, :), uncut_xx_smoothed(firstNonNaNIndex:end), uncut_tvec(firstNonNaNIndex:end), saveDirParticipant, uniqueTable.Patient_ID{curIdx}, init_params, 1, 0)
        % else
        %     plot_bifurcation(dd_orig, uncut_xx_smoothed, uncut_tvec, saveDirParticipant, uniqueTable.Patient_ID{curIdx}, init_params, 1, 0)
        % end

        % Adjust the bifurcation
        [params_tuned_coarsely,rsq_init1,rsq_final1,dd1,xini,iffail1] = tunebif_param([], xx_smoothed, tvec);

        % Further finetune it
        [params_tuned,rsq_ini,mse_ini, rsq_final,mse_final,dd,iffail] = finetune_bif_rev(params_tuned_coarsely,xx_smoothed,tvec,xini, 250);
        
        % Plot the adjusted bifurcation
        [criticalS, c3sol] = plot_bifurcation(dd, xx_smoothed, tvec, saveDirParticipant, participantName, params_tuned, 0, 1);

        if isempty(criticalS)
            continue
        end

        criticalSValue = max(criticalS);
        all_current_hypnogram = participantSubset.CleanSleepStagesInt{curIdx};
        sleepOnsetIdx = participantSubset.SleepOnsetIndex(curIdx);
        time2sleepOnset = participantSubset.Time2Sleep(curIdx);
        hypnogram_cut_post_sleep = all_current_hypnogram(1:sleepOnsetIdx+19);
        if length(hypnogram_cut_post_sleep)>80
            hypnogram_cut = hypnogram_cut_post_sleep(length(hypnogram_cut_post_sleep) - 79: end);
        else 
            hypnogram_cut = hypnogram_cut_post_sleep;
        end 


        [finalCrossing, realCrossingValues, realCrossingTimePoints, crossingValues, crossingTimePoints, crossingDetails, lastPositiveCrossing] = identify_crossing_points(xx_smoothed, criticalSValue, tvec, hypnogram_cut, saveDirParticipant, 0, 'NightNumber',  curIdx);
        
        falseAlarmRate = NaN;
        
        proportionFalseAlarmN1N2 = NaN;
        if ~isempty(realCrossingTimePoints)
            falseAlarmRate = length(realCrossingTimePoints) - 1;
            false_alarm_mask = and([crossingDetails.isRealCrossing], (~[ crossingDetails.isFinalPrediction]));
            ifN1N2 = [crossingDetails.ifN1N2];
            ifFalseAlarmN1N2 = ifN1N2(false_alarm_mask); 
            proportionFalseAlarmN1N2 = sum(ifFalseAlarmN1N2) / length(ifFalseAlarmN1N2);   
            
        end
        %% Log in all the results
        fittingStatsTableRow.Patient_ID = participantName;
        fittingStatsTableRow.Max_Number_Nights_Available(1) = length(non_empty_nights);
        fittingStatsTableRow.Night_Number(1) =  curIdx;
        fittingStatsTableRow.Local_S_Variable{1} = xx_smoothed;
        fittingStatsTableRow.Hypnogram{1} = hypnogram_cut;
        fittingStatsTableRow.UpdatedTimeVec{1} = tvec;
        fittingStatsTableRow.R_square(1) = rsq_final;
        fittingStatsTableRow.MSE(1) = mse_final;
        fittingStatsTableRow.Optimal_Params{1} = params_tuned;
        fittingStatsTableRow.Fitted_Bifurcation{1} = dd(:, 1);
        fittingStatsTableRow.Control_Parameter{1} = dd(:, 2);
        fittingStatsTableRow.SVariableCritical(1) = criticalSValue;
        fittingStatsTableRow.C_3_Solutions{1} = c3sol;
        fittingStatsTableRow.ifNoCrossingPoints(1) = isempty(realCrossingTimePoints);
        fittingStatsTableRow.FinalCrossingTime(1) = finalCrossing; 
        fittingStatsTableRow.CrossingDetails{1} = crossingDetails;
        fittingStatsTableRow.FalseAlarmRate(1) = falseAlarmRate;
        fittingStatsTableRow.FalseAlarmN1N2Proportion(1) = proportionFalseAlarmN1N2;

     fittingStatsTable = [fittingStatsTable; fittingStatsTableRow];
        
    end
   
   save(filePath, 'fittingStatsTable');

end

save(filePath, 'fittingStatsTable');

disp(['Mean of R^2 of fit is ', num2str(mean(fittingStatsTable.R_square, 'omitmissing'))])
disp(['STD of R^2 of fit is ', num2str(std(fittingStatsTable.R_square, 'omitmissing'))])

disp(['Mean of crossing point is ', num2str(mean(fittingStatsTable.FinalCrossingTime, 'omitmissing'))])
disp(['STD of crossing point is ', num2str(std(fittingStatsTable.FinalCrossingTime, 'omitmissing'))])

