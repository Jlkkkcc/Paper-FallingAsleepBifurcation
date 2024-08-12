

num_nights_for_modelling = 1;
saveDir = '/Users/anastasia/Dropbox/ISN work/Work on new dataset/Critical Point Prediction Corrected/';
saveDirAdjustedBifurcation = [saveDir, 'Prediction/', num2str(num_nights_for_modelling), '_Nights_for_Modelling/'];



%% Histogram of prediction times  of the model nights 
handle = figure('Visible', 'off');
histogram(uniqueTable.ModelCrossingTime) %(uniqueTable.ModelCrossingTime <0))

%histogram(earliestTimePoints, 'BinWidth', 0.25); % You can adjust the BinWidth or other properties as needed

% Add a vertical line for the mean
hold on;
xline(mean(uniqueTable.ModelCrossingTime, 'omitnan'), 'r--', 'LineWidth', 2);


%xline( mean(Updated_R_square, 'omitnan'), 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', mean(Updated_R_square, 'omitnan')), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
hold off;

% Add labels and title
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
set(gca, 'Box', 'off')
set(gca, 'FontName', 'Helvetica')
xlabel('Time (min)');
ylabel('Frequency');
%xlim([-10,3])
title(sprintf('Histogram of transition tipping points across training nigths \nMean = %.2f, Std = %.2f',mean(uniqueTable.ModelCrossingTime, 'omitnan'), std(uniqueTable.ModelCrossingTime, 'omitnan')));
grid off;

% Construct filenames for saving
jpegFilename = fullfile(saveDirAdjustedBifurcation, 'Model_Crossing_Times.jpg');
figFilename = fullfile(saveDirAdjustedBifurcation, 'Model_Crossing_Times.fig');
svgFilename = fullfile(saveDirAdjustedBifurcation, 'Model_Crossing_Times.svg');

% Save the plot as JPEG
saveas(handle, jpegFilename);

% Save the figure as .fig file
savefig(handle, figFilename);

% Save the plot as SVG
saveas(handle, svgFilename);

close(handle); % Close the figure



%% Histogram of the R^2 of the model for prediction
handle = figure('Visible', 'off');
histogram(uniqueTable.Updated_R_square)

%histogram(earliestTimePoints, 'BinWidth', 0.25); % You can adjust the BinWidth or other properties as needed

% Add a vertical line for the mean
hold on;
xline(mean(uniqueTable.Updated_R_square, 'omitnan'), 'r--', 'LineWidth', 2);


%xline( mean(Updated_R_square, 'omitnan'), 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', mean(Updated_R_square, 'omitnan')), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
hold off;

% Add labels and title
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
set(gca, 'Box', 'off')
set(gca, 'FontName', 'Helvetica')
xlabel('R^2');
ylabel('Frequency');
%xlim([-10,3])
title(sprintf('Histogram of R^2 of the predictor model \nMean = %.2f, Std = %.2f',mean(uniqueTable.Updated_R_square, 'omitnan'), std(uniqueTable.Updated_R_square, 'omitnan')));
grid off;

% Construct filenames for saving
jpegFilename = fullfile(saveDirAdjustedBifurcation, 'R_squared_of_predictors.jpg');
figFilename = fullfile(saveDirAdjustedBifurcation, 'R_squared_of_predictors.fig');
svgFilename = fullfile(saveDirAdjustedBifurcation, 'R_squared_of_predictors.svg');

% Save the plot as JPEG
saveas(handle, jpegFilename);

% Save the figure as .fig file
savefig(handle, figFilename);

% Save the plot as SVG
saveas(handle, svgFilename);

close(handle); % Close the figure


%% Histogram of the False Alarm Rate
handle = figure('Visible', 'off');
false_alarms =[];
for i = 1:height(uniqueTable)
    false_alarms = [false_alarms, uniqueTable.AllFalseAlarmRates{i}];
end 
histogram(false_alarms)

%histogram(earliestTimePoints, 'BinWidth', 0.25); % You can adjust the BinWidth or other properties as needed

% Add a vertical line for the mean
hold on;
xline(mean(false_alarms, 'omitnan'), 'r--', 'LineWidth', 2);


%xline( mean(Updated_R_square, 'omitnan'), 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', mean(Updated_R_square, 'omitnan')), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
hold off;

% Add labels and title
xlabel('False Alarms per night');
ylabel('Frequency');
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
set(gca, 'Box', 'off')
set(gca, 'FontName', 'Helvetica')
%xlim([-10,3])
title(sprintf('Histogram of False Alarm per predicted night \nMean = %.2f, Std = %.2f',mean(false_alarms, 'omitnan'), std(false_alarms, 'omitnan')));
grid off;


% Construct filenames for saving
jpegFilename = fullfile(saveDirAdjustedBifurcation, 'FalseAlarmRates.jpg');
figFilename = fullfile(saveDirAdjustedBifurcation, 'FalseAlarmRate.fig');
svgFilename = fullfile(saveDirAdjustedBifurcation, 'FalseAlarmRates.svg');

% Save the plot as JPEG
saveas(handle, jpegFilename);

% Save the figure as .fig file
savefig(handle, figFilename);

% Save the plot as SVG
saveas(handle, svgFilename);

close(handle); % Close the figure

%% Histogram of the Mean False Alarm Rate per participant
handle = figure('Visible', 'off');
patients_list = unique(uniqueTable.Patient_ID);
meanFalseAlarmParticipants = [];
for p = 1:length(patients_list)
    patient = patients_list{p};
    patient_table = uniqueTable(strcmp(uniqueTable.Patient_ID, patient), :);
    FalseAlarmParticipant = mean(patient_table.MeanFalseAlarmRate, 'omitmissing');
    meanFalseAlarmParticipants = [meanFalseAlarmParticipants, FalseAlarmParticipant];
end 

histogram(meanFalseAlarmParticipants, 'BinWidth', 0.25);
%histogram(earliestTimePoints, 'BinWidth', 0.25); % You can adjust the BinWidth or other properties as needed

% Add a vertical line for the mean
hold on;
xline(mean(meanFalseAlarmParticipants, 'omitnan'), 'r--', 'LineWidth', 2);


%xline( mean(Updated_R_square, 'omitnan'), 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', mean(Updated_R_square, 'omitnan')), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
hold off;

% Add labels and title
xlabel('Mean False Alarms per participant');
ylabel('Frequency');
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
set(gca, 'Box', 'off')
set(gca, 'FontName', 'Helvetica')
%xlim([-10,3])
title(sprintf('Mean False Alarm per Participant \nMean = %.2f, Std = %.2f',mean(meanFalseAlarmParticipants, 'omitnan'), std(meanFalseAlarmParticipants, 'omitnan')));
grid off;


% Construct filenames for saving
jpegFilename = fullfile(saveDirAdjustedBifurcation, 'MeanParticipantFalseAlarmRates.jpg');
figFilename = fullfile(saveDirAdjustedBifurcation, 'MeanParticipantFalseAlarmRate.fig');
svgFilename = fullfile(saveDirAdjustedBifurcation, 'MeanParticipantFalseAlarmRates.svg');

% Save the plot as JPEG
saveas(handle, jpegFilename);

% Save the figure as .fig file
savefig(handle, figFilename);

% Save the plot as SVG
saveas(handle, svgFilename);

close(handle); % Close the figure


%% Histogram of the SO detection time in predictor model.
handle = figure('Visible', 'off');
patients_list = unique(uniqueTable.Patient_ID);
meanFalseAlarmParticipants = [];
for p = 1:length(patients_list)
    patient = patients_list{p};
    patient_table = uniqueTable(strcmp(uniqueTable.Patient_ID, patient), :);
    FalseAlarmParticipant = mean(patient_table.MeanFalseAlarmRate, 'omitmissing');
    meanFalseAlarmParticipants = [meanFalseAlarmParticipants, FalseAlarmParticipant];
end 

histogram(meanFalseAlarmParticipants, 'BinWidth', 0.25);
%histogram(earliestTimePoints, 'BinWidth', 0.25); % You can adjust the BinWidth or other properties as needed

% Add a vertical line for the mean
hold on;
xline(mean(meanFalseAlarmParticipants, 'omitnan'), 'r--', 'LineWidth', 2);


%xline( mean(Updated_R_square, 'omitnan'), 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', mean(Updated_R_square, 'omitnan')), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'middle');
hold off;

% Add labels and title
xlabel('Mean False Alarms per participant');
ylabel('Frequency');
set(gca,'FontSize', 16)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
set(gca, 'Box', 'off')
set(gca, 'FontName', 'Helvetica')
%xlim([-10,3])
title(sprintf('Mean False Alarm per Participant \nMean = %.2f, Std = %.2f',mean(meanFalseAlarmParticipants, 'omitnan'), std(meanFalseAlarmParticipants, 'omitnan')));
grid off;


% Construct filenames for saving
jpegFilename = fullfile(saveDirAdjustedBifurcation, 'MeanParticipantFalseAlarmRates.jpg');
figFilename = fullfile(saveDirAdjustedBifurcation, 'MeanParticipantFalseAlarmRate.fig');
svgFilename = fullfile(saveDirAdjustedBifurcation, 'MeanParticipantFalseAlarmRates.svg');

% Save the plot as JPEG
saveas(handle, jpegFilename);

% Save the figure as .fig file
savefig(handle, figFilename);

% Save the plot as SVG
saveas(handle, svgFilename);

close(handle); % Close the figure

%% See which variables are the most correlated

% Initialize a structure to hold the correlation coefficients
correlationResults = struct();

% Get the list of variable names
varNames = uniqueTable.Properties.VariableNames;

varNames2Omit = {'MeanFalseAlarmRate', 'STDFalseAlarmRate', 'Mean_EarliestPredictedTimePointLatencies', 'Mean_Predicted_Critical_TimePoints', 'Mean_AverageCriticalityLatency', ...
    'STD_Predicted_Critical_TimePoints', 'STD_EarliestPredictedTimePointLatencies', 'STD_AverageCriticalityLatency'};


% Initialize arrays to hold the results
correlationCoefficients = [];
pValues = [];
significantVariables = {};
% Calculate the correlation coefficient for each variable with MeanFalseAlarmRate
for i = 1:length(varNames)
    varName = varNames{i};
   
    if any(strcmp(varName, varNames2Omit))
        continue
    end 
    % Check if the variable is numeric
    if isnumeric(uniqueTable.(varName)) && ~all(isnan(uniqueTable.(varName)))
        % Calculate the correlation coefficient
        % Calculate the correlation coefficient and p-value
        [R, P] = corr(uniqueTable.(varName), uniqueTable.MeanFalseAlarmRate, 'Type', 'Pearson');
        if isnan(R)
            continue
        end
   
       % Store the results
        correlationCoefficients = [correlationCoefficients; R];
        pValues = [pValues; P];
        significantVariables = [significantVariables; varName];

    end
end

% Create a table to store the results
resultsTable = table(significantVariables, correlationCoefficients, pValues, ...
                     'VariableNames', {'Variable', 'CorrelationCoefficient', 'pValue'});

% Sort the table by the absolute value of the correlation coefficient in descending order
sortedResultsTable = sortrows(resultsTable, 'CorrelationCoefficient', 'descend');

% Display the sorted results table
disp('Variables sorted by their correlation with MeanFalseAlarmRate:');
disp(sortedResultsTable);

% Identify significant correlations (e.g., p-value < 0.05)
significantResults = sortedResultsTable(sortedResultsTable.pValue < 0.05, :);

% Sort significant correlations in ascending order of correlation coefficients
significantResults = sortrows(significantResults, 'CorrelationCoefficient', 'ascend');

% Display significant correlations
disp('Significant correlations (p-value < 0.05):');
disp(significantResults);

% Save the sorted results table to a file
writetable(sortedResultsTable, [saveDirAdjustedBifurcation, 'sortedCorrelationTable.csv']);

% Save the significant results table to a file
writetable(significantResults,  [saveDirAdjustedBifurcation, 'significantCorrelationTable.csv']);





% Initialize arrays to hold the results
correlationCoefficients = [];
pValues = [];
significantVariables = {};
varNames = uniqueTable.Properties.VariableNames;

% Concatenate cell arrays into single numeric arrays

% Handle the concatenation of cell arrays
if iscell(uniqueTable.AllFalseAlarmRates)
    allFalseAlarmRate = [];
    for i = 1:length(uniqueTable.AllFalseAlarmRates)
        if isnumeric(uniqueTable.AllFalseAlarmRates{i})
            allFalseAlarmRate = [allFalseAlarmRate; uniqueTable.AllFalseAlarmRates{i}(:)];
        end
    end
else
    allFalseAlarmRate = uniqueTable.AllFalseAlarmRates;
end

for i = 1:length(varNames)
    varName = varNames{i};
    
    % Check if the variable is numeric or a cell array of numeric arrays
    if isnumeric(uniqueTable.(varName)) || strcmp(varName, 'AllFalseAlarmRates')
        continue
    elseif iscell(uniqueTable.(varName)) && all(cellfun(@isnumeric, uniqueTable.(varName)))
        % Concatenate cell array contents
        variableData = [];
        for j = 1:length(uniqueTable.(varName))
            if isempty(uniqueTable.Updated_C_3_Solutions{j})
                continue
            end 
            variableData = [variableData; uniqueTable.(varName){j}(:)];
        end
    else
        % Skip non-numeric and non-cell array variables
        continue;
    end

    % Ensure the lengths of the concatenated data and allFalseAlarmRate match
    if length(variableData) == length(allFalseAlarmRate)
        % Remove NaNs
        validIdx = ~isnan(variableData) & ~isnan(allFalseAlarmRate);
        variableData = variableData(validIdx);
        allFalseAlarmRateValid = allFalseAlarmRate(validIdx);

        % Calculate the correlation coefficient and p-value
        [R, P] = corr(variableData, allFalseAlarmRateValid, 'Type', 'Pearson');
        
        % Store the results only if the correlation is not NaN
        if ~isnan(R)
            correlationCoefficients = [correlationCoefficients; R];
            pValues = [pValues; P];
            significantVariables = [significantVariables; varName];
        end
    end
end

% Create a table to store the results
resultsTable = table(significantVariables, correlationCoefficients, pValues, ...
                     'VariableNames', {'Variable', 'CorrelationCoefficient', 'pValue'});

% Sort the table by the correlation coefficient in ascending order
sortedResultsTable = sortrows(resultsTable, 'CorrelationCoefficient', 'ascend');

% Display the sorted results table
disp('Variables sorted by their correlation with AllFalseAlarmRate:');
disp(sortedResultsTable);

% Identify significant correlations (e.g., p-value < 0.05)
significantResults = sortedResultsTable(sortedResultsTable.pValue < 0.05, :);

% Display significant correlations
disp('Significant correlations (p-value < 0.05):');
disp(significantResults);

% Save the sorted results table to a file
writetable(sortedResultsTable, [saveDirAdjustedBifurcation, 'sortedCorrelationTable.csv']);

% Save the significant results table to a file
writetable(significantResults,  [saveDirAdjustedBifurcation, 'significantCorrelationTable.csv']);


% Create a bar plot for variables with significant correlations
figure;
bar(categorical(significantResults.Variable, significantResults.Variable), significantResults.CorrelationCoefficient);
ylabel('Correlation Coefficient');
title('Significant Correlations with AllFalseAlarmRate');
xtickangle(45); % Rotate x-axis labels for better readability

% Add p-value annotations to the bar plot
hold on;
for i = 1:height(significantResults)
    text(i, significantResults.CorrelationCoefficient(i), sprintf('p=%.2f', significantResults.pValue(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
hold off;

% Save the bar plot as a figure
saveas(gcf, [saveDirAdjustedBifurcation, 'significantCorrelationsBarPlot.jpg']);
savefig(gcf, [saveDirAdjustedBifurcation, 'significantCorrelationsBarPlot.fig']);
saveas(gcf, [saveDirAdjustedBifurcation, 'significantCorrelationsBarPlot.svg']);

% Create correlation plots for each significant variable with AllFalseAlarmRate
for i = 1:height(significantResults)
    varName = significantResults.Variable{i};
    
    % Extract the data for the variable
    if iscell(uniqueTable.(varName))
        variableData = [];
        for j = 1:length(uniqueTable.(varName))
            if isempty(uniqueTable.Updated_C_3_Solutions{j})
                continue
            end 
               
            variableData = [variableData; uniqueTable.(varName){j}(:)];
             
        end
    else
        variableData = uniqueTable.(varName);
    end
    
    validIdx = ~isnan(variableData) & ~isnan(allFalseAlarmRate);
    variableData = variableData(validIdx);
    allFalseAlarmRateValid = allFalseAlarmRate(validIdx);
    % Create the scatter plot
    figure;
    scatter(variableData, allFalseAlarmRateValid, 'filled');
    hold on;
    
    % Plot the linear fit
    coeffs = polyfit(variableData, allFalseAlarmRateValid, 1);
    fittedX = linspace(min(variableData), max(variableData), 200);
    fittedY = polyval(coeffs, fittedX);
    plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
    
    % Add labels and title
    xlabel(varName);
    ylabel('AllFalseAlarmRate');
    title(['Correlation between ', varName, ' and AllFalseAlarmRate']);
    legend('Data points', 'Linear fit');
    ylim([-0.5, 4.5])
    % Add labels and title
    set(gca,'FontSize', 16)
    set(gca,'TickDir','out')
    set(gca,'ticklength',2*get(gca,'ticklength'))
    set(gca,'lineWidth',2)
    set(gca, 'Box', 'off')
    % Save the plot as a figure
    saveas(gcf, [saveDirAdjustedBifurcation, 'correlation_', varName, '_AllFalseAlarmRate.jpg']);
    saveas(gcf, [saveDirAdjustedBifurcation, 'correlation_', varName, '_AllFalseAlarmRate.svg']);
    savefig(gcf, [saveDirAdjustedBifurcation, 'correlation_', varName, '_AllFalseAlarmRate.fig']);
    close(gcf);
end



% Choose the two variables you want to test against AllFalseAlarmRate
varName1 = 'FluctuationIndex'; % Replace with actual variable name
varName2 = 'RMSSD_S'; % Replace with actual variable name

% Extract the data for the variables
if iscell(uniqueTable.(varName1))
    data1 = [];
    for j = 1:length(uniqueTable.(varName1))
        if isempty(uniqueTable.Updated_C_3_Solutions{j})
            continue
        end 
        data1 = [data1; uniqueTable.(varName1){j}(:)];
    end
else
    data1 = uniqueTable.(varName1);
end

if iscell(uniqueTable.(varName2))
    data2 = [];
    for j = 1:length(uniqueTable.(varName2))
        if isempty(uniqueTable.Updated_C_3_Solutions{j})
            continue
        end 
        data2 = [data2; uniqueTable.(varName2){j}(:)];
    end
else
    data2 = uniqueTable.(varName2);
end

% Ensure that the lengths match
%allFalseAlarmRate = cell2mat(uniqueTable.AllFalseAlarmRates);
validIdx = ~isnan(data1) & ~isnan(data2) & ~isnan(allFalseAlarmRate);
data1 = data1(validIdx);
data2 = data2(validIdx);
allFalseAlarmRate = (allFalseAlarmRate(validIdx'))';

% Create table with appropriate variable names
tbl = table(data1, data2, allFalseAlarmRate', 'VariableNames', {varName1, varName2, 'AllFalseAlarmRate'});

% Perform multiple linear regression
mdl = fitlm(tbl, 'AllFalseAlarmRate ~ FluctuationIndex + RMSSD_S');

% Display model summary
disp(mdl);

% Extract coefficients and p-values
coefficients = mdl.Coefficients.Estimate;
pValues = mdl.Coefficients.pValue;

% Print coefficients and p-values
fprintf('Coefficients:\n');
disp(coefficients);
fprintf('P-values:\n');
disp(pValues);

% Create a 3D scatter plot
figure;
scatter3(data1, data2, allFalseAlarmRate, 'filled');
hold on;

% Generate grid for the fitted plane
[data1Grid, data2Grid] = meshgrid(linspace(min(data1), max(data1), 50), linspace(min(data2), max(data2), 50));
fittedPlane = coefficients(1) + coefficients(2) * data1Grid + coefficients(3) * data2Grid;

% Plot the fitted plane
surf(data1Grid, data2Grid, fittedPlane, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
xlabel(varName1);
ylabel(varName2);
zlabel('AllFalseAlarmRate');
title(['3D Correlation between ', varName1, ', ', varName2, ' and AllFalseAlarmRate']);
legend('Data points', 'Fitted plane');
hold off;

% Save the 3D plot as a figure
saveas(gcf, [saveDirAdjustedBifurcation, '3D_correlation_', varName1, '_', varName2, '_AllFalseAlarmRate.jpg']);
savefig(gcf, [saveDirAdjustedBifurcation, '3D_correlation_', varName1, '_', varName2, '_AllFalseAlarmRate.fig']);