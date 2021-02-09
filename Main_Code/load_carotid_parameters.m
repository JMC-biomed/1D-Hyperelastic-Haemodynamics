
%% Initialize variables.
filename = 'parameter_range.dat';
delimiter = ' ';

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[3,4,5,6,7,8,9,10,11]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [3,4,5,6,7,8,9,10,11]);
rawStringColumns = string(raw(:, [1,2]));


%% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
idx = (rawStringColumns(:, 2) == "<undefined>");
rawStringColumns(idx, 2) = "";

%% Create output variable
parameterrange = table;
parameterrange.I = rawStringColumns(:, 1);
parameterrange.l = categorical(rawStringColumns(:, 2));
parameterrange.c = cell2mat(rawNumericColumns(:, 1));
parameterrange.k1 = cell2mat(rawNumericColumns(:, 2));
parameterrange.k2 = cell2mat(rawNumericColumns(:, 3));
parameterrange.beta = cell2mat(rawNumericColumns(:, 4));
parameterrange.alpha = cell2mat(rawNumericColumns(:, 5));
parameterrange.VarName8 = cell2mat(rawNumericColumns(:, 6));
parameterrange.VarName9 = cell2mat(rawNumericColumns(:, 7));
parameterrange.VarName10 = cell2mat(rawNumericColumns(:, 8));
parameterrange.VarName11 = cell2mat(rawNumericColumns(:, 9));

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns idx;