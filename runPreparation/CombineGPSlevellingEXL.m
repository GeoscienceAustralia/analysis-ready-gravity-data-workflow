clear all
% Folder path
folder = 'Data/GPS_LEVELLING/AUSGeoid2020 Jurisdiction Supplied Data';

% Get list of all Excel files
files = dir(fullfile(folder, '*.xls*'));  % handles .xls and .xlsx

% Initialize empty table to store all data
allData = table();

% Loop through each file
for i = 1:length(files)
    
    % Full file path
    filename = fullfile(folder, files(i).name);
    
    % Read Excel file into table
    T = readtable(filename,'VariableNamingRule','preserve');
    
    % --- Display variable names (first file only, for debugging) ---
    if i == 1
        disp('Column names:');
        disp(filename);
        disp(T.Properties.VariableNames')
    end
    
    state=T.("Jurisdiction");
    role=T.("Nominated Role (modelling, testing)");
    lat = T.("Latitude (DD.DDD)");
    lon = T.("Longitude (DDD.DDD)");
    hGPS= T.("Ellipsoidal Height");                               
    hGPSError= T.("Ellipsoidal height standard uncertainty (1 sigma)");
    ahd= T.("AHD height");
    ahdError  = T.("AHD standard uncertainty (1 sigma)");
  
    if iscell(lon)
        lon = str2double(lon);
    else
        lon = double(lon);
    end

    %Create standardized table
    data = table(state,role,lat, lon, hGPS,hGPSError, ahd, ahdError, ...
        'VariableNames', {'state','role','lat', 'lon', 'hGPS','hGPSError', 'ahd', 'ahdError'});
    
    % Append to master table
    allData = [allData; data];
    
end

disp(allData(1:10,:))

GPS = [allData.lon, ...
       allData.lat, ...
       allData.hGPS - allData.ahd];

% Find rows where h-H is NaN / non-value
idxBad = isnan(GPS(:,3));

% Count them
nBad = sum(idxBad);

% Report
fprintf('Number of rows with non-values in h-H: %d\n', nBad);

% Remove them
GPS_clean = GPS(~idxBad,:);

% Report remaining rows
fprintf('Number of rows remaining after removal: %d\n', size(GPS_clean,1));

% Save clean 3-column matrix to MAT-file
save('Data/GPS_LEVELLING/GPSLevelling2026july.mat','GPS_clean')

 
% Statistics
nPts = size(GPS_clean,1);
meanN = mean(GPS_clean(:,3));
stdN  = std(GPS_clean(:,3));
minN  = min(GPS_clean(:,3));
maxN  = max(GPS_clean(:,3));

figure('Name','GPS Levelling Points (h-H)','NumberTitle','off')

worldmap([-45 -9],[112 155])

load coastlines
plotm(coastlat,coastlon,'k')

scatterm(GPS_clean(:,2),GPS_clean(:,1),8,GPS_clean(:,3),'filled')

colormap(turbo)

cb = colorbar;
cb.Label.String = 'h - H (m)';

title('GPS Levelling Points (h-H)')

% Statistics text box
statsText = sprintf([ ...
    'Count: %d\n' ...
    'Mean: %.4f m\n' ...
    'Std Dev: %.4f m\n' ...
    'Min: %.4f m\n' ...
    'Max: %.4f m'], ...
    nPts, meanN, stdN, minN, maxN);

annotation('textbox',[0.12 0.12 0.22 0.15], ...
    'String',statsText, ...
    'FitBoxToText','on', ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'FontSize',10);
saveas(gcf,'GPSlevellingJuly2026.png') 





