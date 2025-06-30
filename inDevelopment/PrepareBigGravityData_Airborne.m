
filename = 'Data/GRAVITY/AIRBORNE/TransferNow-lastfilesGMEV/GRAV.DAT';  % path to your 13GB file

% Open the file
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open the file.');
end

% Display message
disp('Reading file. This may take some time...');

% Build format string: skip unwanted columns using '%*f'
% '%f' for columns 20, 21, 26, 72, '%*f' for others
numCols = 172;
targetCols = [20, 21, 26, 72, 73, 79, 80, 81];

formatSpec = repmat('%*f ', 1, numCols);
for i = 1:length(targetCols)
    formatSpec = replaceColumnFormat(formatSpec, targetCols(i), '%f ');
end

% Read the file
data = textscan(fid, formatSpec, 'Delimiter', '', 'CollectOutput', true);

% Close the file
fclose(fid);

% Extracted data matrix
selectedData = data{1};

% Assign to separate variables (optional)
longitude = selectedData(:,1);  
latitude = selectedData(:,2); 
elev_MSL = selectedData(:,3); 
GrvFAL56s_GEO_V = selectedData(:,4); 
GrvFAL100s_GEO_V = selectedData(:,5);
FA5000_GEO_V = selectedData(:,6);
BG5000_V= selectedData(:,7);
BG5000_GEO_V= selectedData(:,8);


disp('Done.');

figure
scatter(longitude,latitude,1,GrvFAL56s_GEO_V)
colorbar
colormap(jet)
title(colorbar,'um/s2','FontSize',10);
title('Free air gravity (Geoid), 56sline filter, levelled, vertical')

figure
scatter(longitude,latitude,1,GrvFAL100s_GEO_V)
colorbar
colormap(jet)
title(colorbar,'um/s2','FontSize',10);
title('Free air gravity (Geoid), 100sline filter, levelled, vertical')

FA5000_GEO_V(FA5000_GEO_V == -9999999.99) = NaN;

figure
scatter(longitude,latitude,1,FA5000_GEO_V)
colorbar
colormap(jet)
title(colorbar,'um/s2','FontSize',10);
title('Free Air Gravity, 5000 m full-wavelength spatial filter (Geoid), vertical')

BG5000_V(BG5000_V == -9999999.99) = NaN;

figure
scatter(longitude,latitude,1,-BG5000_V)
caxis([-1210.385 563.724])
colorbar
colormap(jet)
title(colorbar,'um/s2','FontSize',10);
title('Topographic Gravity, 5000 m spatial filter, 2670 kg/m3, vertical')

BG5000_GEO_V(BG5000_GEO_V == -9999999.99) = NaN;

figure
scatter(longitude,latitude,1,-BG5000_GEO_V)
colorbar
colormap(jet)
title(colorbar,'um/s2','FontSize',10);
title('Topographic Gravity, 5000 m full-wavelength spatial filter, 2670 kg/m3 (Geoid), vertical')

%% Helper function to replace %*f with %f at desired column
function fmtOut = replaceColumnFormat(fmtIn, colIndex, newFmt)
    tokens = strsplit(fmtIn);
    tokens{colIndex} = newFmt;
    fmtOut = strjoin(tokens, ' ');
end






