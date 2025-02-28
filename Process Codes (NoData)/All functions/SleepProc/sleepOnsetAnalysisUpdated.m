%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Anastasia Ilina

%% Processes patients hypnogram data identify sleep onset and extract sleep stages prior to sleep onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [firstFiveScores, thirtyBeforeOnsetScores, ifcleanonset, time2sleep, onsetIndex, adjustedStartIndex, ifRecordingStartIsShifted, new_scores, isUnconventionalOnset, hasREMbeforeOnset, sleepOnsetPattern] = sleepOnsetAnalysisUpdated(scores, minutes_before_onset, minutes_after_start, epoch_duration)

    % Constants
    epochs_per_minute = 60 / epoch_duration;  % Since each epoch is 30 seconds

    % Clean the scores from NaNs at the start
    firstNonNanIndex = find(~isnan(scores), 1);
    if isempty(firstNonNanIndex)
        error('Scores are all NaNs');
    end
    
    new_scores = scores(firstNonNanIndex:end);
    adjustedStartIndex = firstNonNanIndex;
    ifRecordingStartIsShifted = firstNonNanIndex > 1;

    % Replace NaNs with a temporary unique value (-1) for onset detection
    scoresForOnsetDetection = new_scores;
    scoresForOnsetDetection(isnan(new_scores)) = -1;

    % Define valid sleep onset patterns
    primaryOnsetPattern = [2, 2];  % Always valid
    secondaryOnsetPatterns = {[2, 3], [2, 4], [3, 3], [3, 4], [4, 4]}; % Only if no prior onset

    onsetIndex = [];
    isUnconventionalOnset = false;
    sleepOnsetPattern = [];

    % Check for primary sleep onset pattern [2,2] first
    for i = 1:length(scoresForOnsetDetection) - 1
        if scoresForOnsetDetection(i) == 2 && scoresForOnsetDetection(i+1) == 2
            onsetIndex = i;
            sleepOnsetPattern = primaryOnsetPattern;
            break;
        end
    end

    % If no primary onset found, check secondary patterns
    if isempty(onsetIndex)
        for pattern = secondaryOnsetPatterns
            patternVec = pattern{1}; % Extract pattern
            for i = 1:length(scoresForOnsetDetection) - 1
                if scoresForOnsetDetection(i) == patternVec(1) && scoresForOnsetDetection(i+1) == patternVec(2)
                    onsetIndex = i;
                    sleepOnsetPattern = patternVec;
                    isUnconventionalOnset = true;
                    break;
                end
            end
            if ~isempty(onsetIndex)
                break; % Stop once an onset is found
            end
        end
    end

    % Check if sleep onset was found
    if isempty(onsetIndex)
        error('No sleep onset found in the provided scores');
    end

    % Check if REM (5) occurred before sleep onset
    hasREMbeforeOnset = any(new_scores(1:onsetIndex-1) == 5);

    % Calculate the time to sleep onset in minutes
    time2sleep = onsetIndex / epochs_per_minute;

    % Extract scores for 5 minutes after sleep onset (including onset)
    endIdx = min(length(new_scores), onsetIndex + minutes_after_start * epochs_per_minute - 1);
    firstFiveScores = new_scores(onsetIndex:endIdx);

    % Extract scores from 30 minutes before onset
    startIdx = max(1, onsetIndex - minutes_before_onset * epochs_per_minute);
    thirtyBeforeOnsetScores = new_scores(startIdx:onsetIndex-1);

    % Check for a "clean" onset (no NREM2 before the detected onset)
    ifcleanonset = isempty(find(new_scores(1:onsetIndex-1) == 2, 1));

end

