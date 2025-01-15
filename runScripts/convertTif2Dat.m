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
