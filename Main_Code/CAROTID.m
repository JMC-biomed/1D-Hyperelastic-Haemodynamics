
%% Import the data
[~, ~, raw] = xlsread('Carotid.xls','Sheet1');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,2);
raw = raw(:,[1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
Vessel = data(:,1);
Abbreviation = cellVectors(:,1);
Length = data(:,2);
Daughter1 = data(:,3);
Daughter2 = data(:,4);
Daughter3 = data(:,5);
Beta_start = data(:,6);
Beta_end = data(:,7);
r0_start = data(:,8);
r0_end = data(:,9);
Vnelem = data(:,10);
Terminal = data(:,11);
Fraction = data(:,12);
Valve = data(:,13);
BP = data(:,14);
RA = data(:,15);
Vein = data(:,16);
Type = data(:,17);
Z = data(:,18);
R = data(:,19);
C = data(:,20);

%% Clear temporary variables
clearvars data raw cellVectors;