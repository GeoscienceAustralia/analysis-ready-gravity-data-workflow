
% Define the path to your .gsb file
filename = 'Data/EXISTING_GEOID_MODELS/AGQG_20201120.gsb';

% Read the .gsb file
T = readgeotable(filename);

% Display the contents of the table
disp(T);

% Plot the data
mapshow(T);
title('Geospatial Data from GSB File');










% Specify the path to your .tif file
tifFile = 'Data/EXISTING_GEOID_MODELS/AGQG_20201120.tif';

% Step 1: Load the raster data and spatial reference
[data, R] = readgeoraster(tifFile);

% Display basic information about the data
disp('Spatial Reference Information:');
disp(R); % Display spatial reference metadata
disp(['Data Dimensions: ', num2str(size(data))]); % Display data dimensions

% Step 2: Handle NoData values (optional)
noDataValue = -9999; % Replace this with the actual NoData value if applicable
if any(data(:) == noDataValue)
    data(data == noDataValue) = NaN; % Replace NoData values with NaN
    disp('NoData values replaced with NaN.');
end

% Step 3: Visualize the data
figure;
% Basic visualization using imagesc
subplot(1, 2, 1);
imagesc(data);
colorbar;
title('Basic Visualization of TIF Data');
xlabel('Column Index');
ylabel('Row Index');

% Georeferenced visualization using mapshow
subplot(1, 2, 2);
mapshow(data, R, 'DisplayType', 'surface');
colorbar;
title('Georeferenced Visualization of TIF Data');

% Step 4: Basic analysis
minValue = min(data(:), [], 'omitnan'); % Minimum value, ignoring NaN
maxValue = max(data(:), [], 'omitnan'); % Maximum value, ignoring NaN
meanValue = mean(data(:), 'omitnan');   % Mean value, ignoring NaN
fprintf('Data Statistics:\nMin: %f, Max: %f, Mean: %f\n', minValue, maxValue, meanValue);


% Loop through each band and visualize
numBands = size(data, 3);
figure;
for i = 1:numBands
    subplot(2, 2, i); % Adjust grid size if more bands
    imagesc(data(:, :, i));
    colorbar;
    title(['Band ', num2str(i)]);
    xlabel('Column Index');
    ylabel('Row Index');
end















% Step 1: Specify the folder containing the .tif files
inputFolder = 'C:\analysis-ready-gravity-data-workflow\outputs\GridsvicOld'; % Replace with your folder path
outputFolder = 'C:\analysis-ready-gravity-data-workflow\outputs\GridsvicOld'; % Replace with your desired output folder

% Ensure the output folder exists
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

% Step 2: Get a list of all .tif files in the folder
tifFiles = dir(fullfile(inputFolder, '*.tif'));

% Step 3: Loop through each .tif file
for k = 1:length(tifFiles)
    % Get the full path of the current .tif file
    tifFilePath = fullfile(inputFolder, tifFiles(k).name);
    
    % Read the .tif file
    [tifData, R] = geotiffread(tifFilePath); % For GeoTIFF files
    tifData = double(tifData); % Convert to double precision if necessary
    
    % Get the file name without extension
    [~, fileName, ~] = fileparts(tifFiles(k).name);
    
    % Create the corresponding .DAT file name
    datFilePath = fullfile(outputFolder, [fileName, '.DAT']);
    
    % Step 4: Open the .DAT file for writing
    fileID = fopen(datFilePath, 'w');
    
    % Step 5: Write metadata (if needed)
    fprintf(fileID, 'WINTER .DAT FILE\n');
    fprintf(fileID, 'Rows: %d, Columns: %d\n', size(tifData, 1), size(tifData, 2));
    
    % Step 6: Write data row by row
    for i = 1:size(tifData, 1)
        fprintf(fileID, '%10.3f\t', tifData(i, :)); % Tab-delimited
        fprintf(fileID, '\n'); % Newline for each row
    end
    
    % Close the .DAT file
    fclose(fileID);
    
    % Display progress
    disp(['Converted: ', tifFiles(k).name, ' to ', fileName, '.DAT']);
end

disp('All .tif files have been processed and converted to .DAT format.');
