function [finalCrossing, realCrossingValues, realCrossingTimePoints, crossingValues, crossingTimePoints, crossingDetails, lastPositiveCrossing] = identify_crossing_points(xx_smoothed, criticalSValue, tvec, hypnogram,  save_dir, iffittedmodel, varargin)
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Author: Anastasia Ilina

% Identifies points of the crossing of the S-variable threshold for sleep
% onset prediction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dd = [];
    if ~isempty(varargin)
        for idx = 1:length(varargin)
            if strcmp(varargin{idx},  'NightNumber')
                night_number = varargin{idx+1};
            end 
            if strcmp(varargin{idx}, 'Bifurcation')
                dd = varargin{idx+1};
            end 
            
        end 
    end 
    
    
    starting_length = length(hypnogram)*0.5 - 10;
    % Identify crossing indices
    crossingIndices = find((xx_smoothed(1:end-1) < criticalSValue & xx_smoothed(2:end) >= criticalSValue) | ...
                           (xx_smoothed(1:end-1) > criticalSValue & xx_smoothed(2:end) <= criticalSValue));
    
    
    if hypnogram(length(hypnogram) - 20) == 2
        hypnogram = hypnogram(1:end - 1);
    end 

    
    % Initialize the expanded hypnogram
    expanded_hypnogram = zeros(1, length(xx_smoothed));
    
    % Loop through each hypnogram epoch and broadcast its value
    for i = 1:length(hypnogram)
        % Each hypnogram epoch corresponds to 10 s-variable epochs
        start_idx = (i - 1) * 10 + 1;
        end_idx = start_idx + 10;
        
        % Broadcast the hypnogram value
        expanded_hypnogram(start_idx:end_idx) = hypnogram(i);
    end
     
    if length(expanded_hypnogram) > length(xx_smoothed)
        expanded_hypnogram = expanded_hypnogram(length(expanded_hypnogram) - length(xx_smoothed):end);

    end 

    time2fallasleep = abs((tvec(1)) - 0.1);

    % Initialize output variables
    realCrossingTimePoints = [];
    realCrossingValues = [];
    crossingTimePoints = [];
    crossingValues = [];
    crossingDetails = struct('downwardCrossingTime', [], 'upwardCrossingTime', [], 'timeSpentBelow', [], 'isRealCrossing', [], 'isFinalPrediction', [], 'ifpositive', [], 'sleepStageDownwardCrossing', [], 'sleepStageUpwardCrossing', [], 'sleepStagesDuring',[], 'ifN1N2', []);
    finalCrossing = NaN; 
    lastPositiveCrossing = NaN;


    if isempty(crossingIndices)
        return;
    end
    
    % Interpolate to find the exact crossing time points
    iffinalpredictionmade = false;
    for i = 1:length(crossingIndices)
        idx = crossingIndices(i);
        % Linear interpolation
        t1 = tvec(idx);
        t2 = tvec(idx + 1);
        x1 = xx_smoothed(idx);
        x2 = xx_smoothed(idx + 1);
        % Interpolating to find the crossing time
        timepoint_int = t1 + (criticalSValue - x1) * (t2 - t1) / (x2 - x1);
        ifpositive  = false;
        isFinalPrediction = false;
        % Determine crossing type
        
        if x1 > criticalSValue && x2 <= criticalSValue
            % Downward crossing
            downwardCrossingTime = timepoint_int;
            
            % Find corresponding upward crossing
            upwardIdx = find(tvec(idx+1:end) > downwardCrossingTime & xx_smoothed(idx+1:end) >= criticalSValue, 1);
            
            % Calculate the index in the expanded hypnogram
            downwardIndex = floor((time2fallasleep + downwardCrossingTime) / 0.05);
            
            % Ensure the index is within the bounds of the expanded hypnogram
            downwardIndex = min(max(downwardIndex, 1), length(expanded_hypnogram));
            
            % Identify the sleep stage
            sleep_stage_at_crossing = expanded_hypnogram(downwardIndex);
            
            if downwardCrossingTime > 0
               ifpositive = true;
            end 

            if ~isempty(upwardIdx) 
                upwardCrossingTime = tvec(idx + upwardIdx);
                % Calculate time spent below the threshold
                upwardIndex = floor((time2fallasleep + upwardCrossingTime) / 0.05);
                if upwardIndex > length(expanded_hypnogram)
                    upwardIndex = length(expanded_hypnogram);
                end 
                upwardSleepStage = expanded_hypnogram(upwardIndex);
                hypno_segment = expanded_hypnogram(downwardIndex:upwardIndex);

                % Check if it's a real crossing
                timeSpentBelow = upwardCrossingTime - downwardCrossingTime;
                if upwardCrossingTime > 0 
                    if downwardCrossingTime > 0
                     
                        isRealCrossing = false;
                        
                    else
                        isFinalPrediction = true;
                        iffinalpredictionmade = true;
                        isRealCrossing = true;
                        timeSpentBelow = 0 - downwardCrossingTime;
                    end 
                else
                    isRealCrossing = (timeSpentBelow > 1 || downwardCrossingTime > -1.0);
                end 
            else
                upwardSleepStage = NaN;
                upwardCrossingTime = NaN;
                hypno_segment = [];
                timeSpentBelow = 0 - downwardCrossingTime;
                isRealCrossing = false;
                if ~iffinalpredictionmade
                    isFinalPrediction = true;
                    isRealCrossing = true;
                end 
            end

            
            % Store the details of the crossing
            crossingDetails(end+1) = struct('downwardCrossingTime', downwardCrossingTime, ...
                                            'upwardCrossingTime', upwardCrossingTime, ...
                                            'timeSpentBelow', timeSpentBelow, ...
                                            'isRealCrossing', isRealCrossing,...
                                            'isFinalPrediction', isFinalPrediction,...
                                            'ifpositive', ifpositive, ...
                                            'sleepStageDownwardCrossing', sleep_stage_at_crossing, ...
                                            'sleepStageUpwardCrossing', upwardSleepStage, ...
                                            'sleepStagesDuring',hypno_segment, ...
                                            'ifN1N2', sum(hypno_segment)>0);
        end
    end

    crossingDetails(1) = [];
    
    % Extract crossing time points and values for plotting
    crossingTimePoints = [crossingDetails.downwardCrossingTime];
    crossingValues = repmat(criticalSValue, size(crossingTimePoints));
    
    realCrossingTimePoints = crossingTimePoints(logical([crossingDetails.isRealCrossing]));
    realCrossingValues = repmat(criticalSValue, size(realCrossingTimePoints));
    
    if sum([crossingDetails.ifpositive]) > 0
        lastPositiveCrossing = max(crossingTimePoints(logical([crossingDetails.ifpositive])));
    else
        lastPositiveCrossing = NaN;
    end 
    
    if sum([crossingDetails.isFinalPrediction]) == 0
        finalCrossing = lastPositiveCrossing;
    else
        finalCrossing = crossingTimePoints(logical([crossingDetails.isFinalPrediction]));
    end

    % Plotting the results
    handle = figure('Visible', 'off');
    hold on;
    set(gca,'FontSize', 16);
    set(gca,'TickDir','out');
    set(gca,'TickLength', 2*get(gca,'TickLength'));
    set(gca,'LineWidth', 2);
    set(gca, 'Box', 'off');
    set(gca, 'FontName', 'Helvetica');
    ylabel('s');
    xlabel('Time (min)');
    ax = gca;
    ylim([0, 14]);
    
    plot(tvec, xx_smoothed, 'k-', 'LineWidth', 1.5); % Plot xx_smooth vs tvec
    if iffittedmodel && ~isempty(dd)
        plot(tvec, dd, 'b-', 'LineWidth', 2)
    end 
    yline(criticalSValue, 'g--', 'LineWidth', 1.5); % Plot criticalSValue as a horizontal line
    
    % Mark crossing points
    for i = 1:length(crossingDetails)
        if crossingDetails(i).isRealCrossing
            plot(crossingDetails(i).downwardCrossingTime, criticalSValue, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Plot real crossing points
        else
            plot(crossingDetails(i).downwardCrossingTime, criticalSValue, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Plot other crossing points
        end
    end
    
    xticks(-30:5:10)
    xlim([-30, 10])
    line([0, 0], ax.YLim, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r');
    title('S-variable vs Time with Critical Value and Crossing Points');
    if iffittedmodel
        legend('S-variable', 'Fitted Bifurcation', 'Critical Value', 'Crossing Points');
    else
        legend('S-variable', 'Critical Value', 'Crossing Points');
    end 
    grid off;
    hold off;

    % Construct the title with line breaks
    if ~isempty(crossingTimePoints)
        graphTitle = {['Critical transition detected ', num2str(crossingTimePoints(1)), ' minutes before SO']};
    else
        graphTitle = {'No Critical Transition Detected'};
    end
    
    % Set the title
    title(graphTitle, 'Interpreter', 'none');
    
    % Construct filenames for saving
    if iffittedmodel
        jpegFilename = fullfile(save_dir, 'CriticalTransitionDetected.jpg');
        figFilename = fullfile(save_dir, 'CriticalTransitionDetected.fig');
        svgFilename = fullfile(save_dir, 'CriticalTransitionDetected.svg');
    else
        jpegFilename = fullfile(save_dir, ['Night_', num2str(night_number), 'CriticalTransitionDetected.jpg']);
        figFilename = fullfile(save_dir, ['Night_', num2str(night_number), 'CriticalTransitionDetected.fig']);
        svgFilename = fullfile(save_dir, ['Night_', num2str(night_number), 'CriticalTransitionDetected.svg']);
    end 
    
    % Save the plot as JPEG
    saveas(handle, jpegFilename);
    
    % Save the plot as JPEG
    saveas(handle, svgFilename);
    
    % Save the figure as .fig file
    savefig(handle, figFilename);
    
    close(handle); % Close the figure


    
end
