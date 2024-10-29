%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Anastasia Ilina

%% Processes patients hypnogram data identify sleep onset and extract sleep stages prior to sleep onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [firstFiveScores, thirtyBeforeOnsetScores, ifcleanonset, time2sleep, onsetIndex, adjustedStartIndex, ifRecordingStartIsShifted, new_scores] = sleepOnsetAnalysis(scores, minutes_before_onset, minutes_after_start, epoch_duration)
    
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

    % Identify first occurrence of 2 consecutive epochs scored as "2"
    onsetIndex = find(conv(scoresForOnsetDetection, [1, 1], 'valid') == 4, 1); % '4' because 2+2=4

    % Check if sleep onset was found
    if isempty(onsetIndex)
        error('No sleep onset found in the provided scores');
    end

    % Calculate the time to sleep onset in minutes
    time2sleep = onsetIndex / epochs_per_minute;

    % Extract scores for 5 minutes after sleep onset (including onset)
    endIdx = min(length(new_scores), onsetIndex + minutes_after_start * epochs_per_minute - 1);
    firstFiveScores = new_scores(onsetIndex:endIdx);

    % Extract scores from 30 minutes before onset
    startIdx = max(1, onsetIndex - minutes_before_onset * epochs_per_minute);
    thirtyBeforeOnsetScores = new_scores(startIdx:onsetIndex-1);

    % Check for a "clean" onset
    ifcleanonset = isempty(find(new_scores(1:onsetIndex-1) == 2, 1));

end
