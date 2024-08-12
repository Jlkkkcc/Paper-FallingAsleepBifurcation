function plot_timepoint_histograms(earliestTimePoints, latestTimePoints, meanTimePoints, allTimePoints, savedir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Author: Anastasia Ilina

% Performs ploting of the histograms for different latencies of critical
% point identification

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meanEarliestTimepoint = mean(earliestTimePoints, 'omitmissing');
stdEarliestTimepoint = std(earliestTimePoints, 'omitmissing');
meanLatestTimepoint = mean(latestTimePoints, 'omitmissing');
stdLatestTimepoint = std(latestTimePoints, 'omitmissing');
meanMeanTimepoint = mean(meanTimePoints, 'omitmissing');
stdMeanTimepoint = std(meanTimePoints, 'omitmissing');
meanAllTimepoint = mean(allTimePoints, 'omitmissing');
stdAllTimepoint = std(allTimePoints, 'omitmissing');

%% Histogram of the earliest time points 
handle = figure('Visible', 'off');
histogram(earliestTimePoints, 'BinWidth', 0.25); % You can adjust the BinWidth or other properties as needed

% Add a vertical line for the mean
hold on;
xline(meanEarliestTimepoint, 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', meanEarliestTimepoint), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
hold off;

% Add labels and title
xlabel('Time Points');
ylabel('Frequency');
%xlim([-10,3])
title(sprintf('Histogram of Earliest Time Points for Critical Transition\nMean = %.2f, Std = %.2f', meanEarliestTimepoint, stdEarliestTimepoint));
grid off;

% Construct filenames for saving
jpegFilename = fullfile(savedir, 'EarliestTransitionDetection.jpg');
figFilename = fullfile(savedir, 'EarliestTransitionDetection.fig');

% Save the plot as JPEG
saveas(handle, jpegFilename);

% Save the figure as .fig file
savefig(handle, figFilename);

close(handle); % Close the figure

%% Histogram of the latest time points 
handle = figure('Visible', 'off');
histogram(latestTimePoints, 'BinWidth', 0.25); % You can adjust the BinWidth or other properties as needed

% Add a vertical line for the mean
hold on;
xline(meanLatestTimepoint, 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', meanLatestTimepoint), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');

hold off;

% Add labels and title
xlabel('Time Points');
ylabel('Frequency');
%xlim([-10,3])
title(sprintf('Histogram of Latest Time Points for Critical Transition\nMean = %.2f, Std = %.2f', meanLatestTimepoint, stdLatestTimepoint));
grid off;

% Construct filenames for saving
jpegFilename = fullfile(savedir, 'LatestTransitionDetection.jpg');
figFilename = fullfile(savedir, 'LatestTransitionDetection.fig');

% Save the plot as JPEG
saveas(handle, jpegFilename);

% Save the figure as .fig file
savefig(handle, figFilename);

close(handle); % Close the figure

%% Histogram of the mean of detection timepoints
handle = figure('Visible', 'off');
histogram(meanTimePoints, 'BinWidth', 0.25); % You can adjust the BinWidth or other properties as needed

% Add a vertical line for the mean
hold on;
xline(meanMeanTimepoint, 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', meanMeanTimepoint), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
hold off;

% Add labels and title
xlabel('Time Points');
ylabel('Frequency');
%xlim([-10,3])
title(sprintf('Histogram of Mean Time Points for Critical Transition\nMean = %.2f, Std = %.2f', meanMeanTimepoint, stdMeanTimepoint));
grid off;

% Construct filenames for saving
jpegFilename = fullfile(savedir, 'MeanTransitionDetection.jpg');
figFilename = fullfile(savedir, 'MeanTransitionDetection.fig');

% Save the plot as JPEG
saveas(handle, jpegFilename);

% Save the figure as .fig file
savefig(handle, figFilename);

close(handle); % Close the figure

%% Histogram of all of the timepoints detected 
handle = figure('Visible', 'off');
histogram(allTimePoints, 'BinWidth', 0.25); % You can adjust the BinWidth or other properties as needed

% Add a vertical line for the mean
hold on;
xline(meanAllTimepoint, 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', meanAllTimepoint), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
hold off;

% Add labels and title
xlabel('Time Points');
ylabel('Frequency');
%xlim([-10,3])
title(sprintf('Histogram of All Time Points for Critical Transition\nMean = %.2f, Std = %.2f', meanAllTimepoint, stdAllTimepoint));
grid off;

% Construct filenames for saving
jpegFilename = fullfile(savedir, 'AllTransitionDetection.jpg');
figFilename = fullfile(savedir, 'AllTransitionDetection.fig');

% Save the plot as JPEG
saveas(handle, jpegFilename);

% Save the figure as .fig file
savefig(handle, figFilename);

close(handle); % Close the figure