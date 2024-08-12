function plot_predictions_together(estimatedSVariables, criticalSValue, globalMin, saveDirParticipant)

smoothingMethod = 'retrospective_median';
smoothingWindow = 20;
smoothingOrder = 4;

handle = figure('Visible', 'off');
set(handle, 'Position', [100, 100, 800, 1200]); % [left, bottom, width, height]

% Initialize arrays to store plot handles for the legend
plotHandles = [];
legendLabels = {'S-variable', 'Critical Value', 'Crossing Points'};

% Number of subplots
numSubplots = length(estimatedSVariables);
numRows = 4;
numCols = 2;

for n_n = 1:numSubplots

    distPerNight1 = estimatedSVariables{n_n};
    uncut_xx_smoothed = smoothdata(distPerNight1, 'movmedian', [smoothingWindow, 0], 'omitmissing');
    uncut_xx_smoothed  = uncut_xx_smoothed - globalMin;
    uncut_xx_smoothed = uncut_xx_smoothed(1:length(uncut_xx_smoothed) - 9);
    uncut_tvec =  (-30*60 + 6):3:10*60;
    uncut_tvec = uncut_tvec / 60;
    tvec_night = uncut_tvec(max(1, length(uncut_tvec) - length(uncut_xx_smoothed)+1): end);
    xx_smoothed = uncut_xx_smoothed(max(1,  length(uncut_xx_smoothed) - length(tvec_night) +1): end);

    % Identify crossing indices
    crossingIndices = find((xx_smoothed(1:end-1) < criticalSValue & xx_smoothed(2:end) >= criticalSValue) | ...
        (xx_smoothed(1:end-1) > criticalSValue & xx_smoothed(2:end) <= criticalSValue));

    % Initialize output variables
    realCrossingTimePoints = [];
    realCrossingValues = [];
    crossingTimePoints = [];
    crossingValues = [];
    crossingDetails = struct('downwardCrossingTime', [], 'upwardCrossingTime', [], 'timeSpentBelow', [], 'isRealCrossing', []);

    if isempty(crossingIndices)
        return;
    end

    % Interpolate to find the exact crossing time points
    for i = 1:length(crossingIndices)
        idx = crossingIndices(i);
        % Linear interpolation
        t1 = tvec_night(idx);
        t2 = tvec_night(idx + 1);
        x1 = xx_smoothed(idx);
        x2 = xx_smoothed(idx + 1);
        % Interpolating to find the crossing time
        timepoint_int = t1 + (criticalSValue - x1) * (t2 - t1) / (x2 - x1);

        if timepoint_int > 0
            continue
        end
        % Determine crossing type
        if x1 > criticalSValue && x2 <= criticalSValue
            % Downward crossing
            downwardCrossingTime = timepoint_int;

            % Find corresponding upward crossing
            upwardIdx = find(tvec_night(idx+1:end) > downwardCrossingTime & xx_smoothed(idx+1:end) >= criticalSValue, 1);
            if ~isempty(upwardIdx)
                upwardCrossingTime = tvec_night(idx + upwardIdx);
                % Calculate time spent below the threshold

                timeSpentBelow = upwardCrossingTime - downwardCrossingTime;
                % Check if it's a real crossing
                isRealCrossing = (timeSpentBelow > 1 || downwardCrossingTime > -1.0);
            else
                upwardCrossingTime = NaN;
                timeSpentBelow = 0 - downwardCrossingTime;
                isRealCrossing = true;
            end

            % Store the details of the crossing
            crossingDetails(end+1) = struct('downwardCrossingTime', downwardCrossingTime, ...
                'upwardCrossingTime', upwardCrossingTime, ...
                'timeSpentBelow', timeSpentBelow, ...
                'isRealCrossing', isRealCrossing);
        end
    end

    crossingDetails(1) = [];

    % Extract crossing time points and values for plotting
    crossingTimePoints = [crossingDetails.downwardCrossingTime];
    crossingValues = repmat(criticalSValue, size(crossingTimePoints));

    realCrossingTimePoints = crossingTimePoints(logical([crossingDetails.isRealCrossing]));
    realCrossingValues = repmat(criticalSValue, size(realCrossingTimePoints));



    % Create a subplot
    subplot(numRows, numCols, n_n);
    hold on;
    %set(gca,'FontSize', 16)
    set(gca,'TickDir','out')
    set(gca,'ticklength',2*get(gca,'ticklength'))
    set(gca,'lineWidth',2)
    set(gca, 'Box', 'off')
    set(gca, 'FontName', 'Helvetica')
    ylabel('s')
    xlim([-30,10])
    xlabel('Time (min)')

    ax = gca;
    ylim([0,14])

    h1 = plot(tvec_night, xx_smoothed, 'k-', 'LineWidth', 1.5); % Plot xx_smooth vs tvec
    h2 = yline(criticalSValue, 'g--', 'LineWidth', 1.5); % Plot criticalSValue as a horizontal line

    % Mark crossing points
    if ~isempty(crossingTimePoints)
    for i = 1:length(crossingDetails)
        if crossingDetails(i).isRealCrossing
            h3 = plot(crossingDetails(i).downwardCrossingTime, criticalSValue, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Plot real crossing points
        else
            h3 = plot(crossingDetails(i).downwardCrossingTime, criticalSValue, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Plot other crossing points
        end
    end
    end 

    % % Mark crossing points
    % if ~isempty(realCrossingTimePoints)
    %     h3 = plot(realCrossingTimePoints, realCrossingValues , 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Plot crossing points
    % end

    line([0,0],ax.YLim,'LineStyle','--','LineWidth',2,'Color','r')
    if n_n == 1
        if ~isempty(crossingTimePoints)
        % Only add to plotHandles if it's the first subplot
            plotHandles = [h1, h2, h3];
        else
             plotHandles = [h1, h2];
        end 
    end
    hold off;

    % Construct the title with line breaks
    if ~isempty( realCrossingTimePoints)
        graphTitle = {['Critical transition detected ', num2str(realCrossingTimePoints(end)), ' minutes before SO']};
    else
        % Construct the title with line breaks

        graphTitle = {'No Critical transitions detected '};
    end

    % Set the title for each subplot
    title(graphTitle, 'Interpreter', 'none');


end

% Add one overall title for the entire figure
sgtitle('S-variable vs Time with Critical Value and Crossing Points for All Nights');

% Create an invisible axis for the legend
legAx = axes('Position', [0, 0, 1, 1], 'Visible', 'off');
%legend(legAx, plotHandles, legendLabels, 'Location', 'southoutside', 'Orientation', 'horizontal');

% Save the subplot figure
jpegFilename = fullfile(saveDirParticipant, 'AllNights_CriticalTransitionDetected.jpg');
figFilename = fullfile(saveDirParticipant, 'AllNights_CriticalTransitionDetected.fig');
svgFilename = fullfile(saveDirParticipant, 'AllNights_CriticalTransitionDetected.svg');

% Save the plot as JPEG
saveas(handle, jpegFilename);

% Save the figure as .fig file
savefig(handle, figFilename);

% Save the plot as JPEG
saveas(handle, svgFilename);

close(handle); % Close the figure