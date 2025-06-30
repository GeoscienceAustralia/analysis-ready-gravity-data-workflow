close all
clear

% Look at Victoria Data.
Vic_Data=importdata('Data\GRAVITY\AIRBORNE/23102024victoriaOtter/FD012_Grav.csv');
Vic_Data.data(end,:)=[];
Vic_LongOtter=round(Vic_Data.data(:,11)*60)/60;
Vic_LatOtter=round(Vic_Data.data(:,10)*60)/60;
Vic_HOtter=Vic_Data.data(:,9);
Vic_Grav_anomOtter=Vic_Data.data(:,36);

figure
scatter(Vic_LongOtter,Vic_LatOtter,1,Vic_Grav_anomOtter)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('Otter')

% Look at Adelaide Data.
% Adelaide_Data=importdata('Data/GRAVITY/AIRBORNE/adelaide2025-06-05/GRAV.DAT');
% Adelaide_Long=round(Adelaide_Data(:,20)*60)/60;
% Adelaide_Lat=round(Adelaide_Data(:,21)*60)/60;
% Adelaide_H=Adelaide_Data(:,26);
% Adelaide_Grav_anom=Adelaide_Data(:,79);
% Adelaide_Grav_anom(Adelaide_Grav_anom == -9999999.99) = NaN;
% Adelaide_Grav_anom=Adelaide_Grav_anom/10;
% 
% figure
% scatter(Adelaide_Long,Adelaide_Lat,1,Adelaide_Grav_anom)
% colorbar
% colormap(jet)
% title(colorbar,'mGal','FontSize',10);
% title('FinalFreeAirGravity,spatialfilter(Geoid),vertical')
% 
% figure
% scatter(Adelaide_Long,Adelaide_Lat,1,Adelaide_H)
% colorbar
% colormap(jet)
% title(colorbar,'m','FontSize',10);
% title('elevMSL')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Look at final Victoria Data.

filename = 'Data/GRAVITY/AIRBORNE/TransferNow-lastfilesGMEV/GRAV.DAT';  % path to your 13GB file

% Open the file
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open the file.');
end

% Display message
disp('Reading file. This may take some time...');

numCols = 171;
targetCols = [20, 21, 26, 73];

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
longitude = selectedData(:,1)*60/60;  
latitude = selectedData(:,2)*60/60; 
elev_MSL = selectedData(:,3); 
FA5000_GEO_V = selectedData(:,4);
FA5000_GEO_V(FA5000_GEO_V == -9999999.99) = NaN;
FA5000_GEO_V=FA5000_GEO_V/10;

disp('Done.');

figure
scatter(longitude,latitude,1, FA5000_GEO_V)
colorbar
colormap(jet)
title(colorbar,'um/s2','FontSize',10);
title('FinalFreeAirGravity,spatialfilter(Geoid),vertical')

figure
scatter(longitude,latitude,1,elev_MSL)
colorbar
colormap(jet)
title(colorbar,'m','FontSize',10);
title('elevMSL')
















% Helper function to replace %*f with %f at desired column
function fmtOut = replaceColumnFormat(fmtIn, colIndex, newFmt)
    tokens = strsplit(fmtIn);
    tokens{colIndex} = newFmt;
    fmtOut = strjoin(tokens, ' ');
end






