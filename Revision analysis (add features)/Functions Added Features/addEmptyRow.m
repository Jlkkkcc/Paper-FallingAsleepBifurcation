function newTable = addEmptyRow(inputTable)
    % Get variable types and names
    varTypes = varfun(@class, inputTable, 'OutputFormat', 'cell');
    varNames = inputTable.Properties.VariableNames;
    
    % Create empty values for each type
    emptyValues = cell(1, width(inputTable));
    
    for i = 1:length(varTypes)
        switch varTypes{i}
            case 'double'
                emptyValues{i} = NaN;
            case 'single'
                emptyValues{i} = single(NaN);
            case 'categorical'
                emptyValues{i} = categorical(missing);
            case 'cell'
                emptyValues{i} = {[]};
            case 'string'
                emptyValues{i} = missing;
            case 'datetime'
                emptyValues{i} = NaT;
            case 'logical'
                emptyValues{i} = false;
            case 'char'
                emptyValues{i} = '';
            case 'struct'
                % Get the structure from the first row of the table
                originalStruct = inputTable{1, i};
                % Create empty struct with same fields
                emptyStruct = structfun(@(x) [], originalStruct, 'UniformOutput', false);
                emptyValues{i} = emptyStruct;
            otherwise
                emptyValues{i} = missing;
        end
    end
    
    % Create empty row
    emptyRow = table(emptyValues{:}, 'VariableNames', varNames);
    
    % Append to original table
    newTable = [inputTable; emptyRow];
end