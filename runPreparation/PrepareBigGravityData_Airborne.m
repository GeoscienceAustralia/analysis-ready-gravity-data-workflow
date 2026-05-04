close all
clear

addpath('functions');

% Look at WA Data.
WA_Data=importdata('Data/GRAVITY/AIRBORNE/WA/WA_Airborne.mat');
WA_Long=round(WA_Data(:,1)*60)/60;
WA_Lat=round(WA_Data(:,2)*60)/60;
WA_H=WA_Data(:,3);
WA_Grav_anom=WA_Data(:,4);

% Look at Adelaide Data.
Adelaide_Data=importdata('Data/GRAVITY/AIRBORNE/adelaide2025-06-05/GRAV.DAT');
Adelaide_Long=round(Adelaide_Data(:,20)*60)/60;
Adelaide_Lat=round(Adelaide_Data(:,21)*60)/60;
Adelaide_H=Adelaide_Data(:,26);
Adelaide_Grav_anom=Adelaide_Data(:,73)/10;

% Look at Gippsland Data.
AirborneGippsland = importdata('Data/GRAVITY/AIRBORNE/Gippsland/GRAV.DAT');
Gippsland_Long=round(AirborneGippsland(:,10)*60)/60;
Gippsland_Lat=round(AirborneGippsland(:,11)*60)/60;
Gippsland_H=AirborneGippsland(:,13);
Gippsland_Grav_anom=AirborneGippsland(:,36);

% Look at final Victoria Data.

filename = 'Data\GRAVITY\AIRBORNE\TransferNow-lastfilesfromGMEV\GRAV.DAT';  % path to your 13GB file

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
GrvFAL100s_GEO_V = selectedData(:,4)/10;

disp('Done.');

% Look at NSW data.
NSW_Data=importdata('Data\GRAVITY\AIRBORNE\18042024victoriaNSW/gravNSW.csv');
NSW_Data.data(end,:)=[];
NSW_Long=round(NSW_Data.data(:,65)*60)/60;
NSW_Lat=round(NSW_Data.data(:,62)*60)/60;
NSW_H=NSW_Data.data(:,66);
NSW_Grav_anom=NSW_Data.data(:,36)/10;% its in mico m's per secon ^2 for some reason. i.e. needs dividing by 10 to get into mGal.

% Check by comparing to EGM2008
GGM=importdata('Data/GGM/EGM2008_For_Gridded_Int.mat');
GGM_Gi=griddedInterpolant(GGM.x,GGM.y,GGM.z,GGM.g);

GGM_Gi_interpolatedVIC=GGM_Gi(longitude,-latitude,elev_MSL);
GGM_Gi_interpolatedAdelaide=GGM_Gi(Adelaide_Long,-Adelaide_Lat,Adelaide_H);
GGM_Gi_interpolatedWA=GGM_Gi(WA_Long,-WA_Lat,WA_H);
GGM_Gi_interpolatedGippsland=GGM_Gi(Gippsland_Long,-Gippsland_Lat,Gippsland_H);
GGM_Gi_interpolatedNSW=GGM_Gi(NSW_Long,-NSW_Lat,NSW_H);

% plots
plotGGMGravityComparison( ...
    WA_Long, WA_Lat, WA_Grav_anom, GGM_Gi_interpolatedWA, 'outputs/plots/EGM2008WA.png');

plotGGMGravityComparison( ...
    longitude, latitude, GrvFAL100s_GEO_V, GGM_Gi_interpolatedVIC, 'outputs/plots/EGM2008Vic.png');

plotGGMGravityComparison( ...
    Adelaide_Long, Adelaide_Lat, Adelaide_Grav_anom, GGM_Gi_interpolatedAdelaide, 'outputs/plots/EGM2008Adelaide.png');

plotGGMGravityComparison( ...
    Gippsland_Long, Gippsland_Lat, Gippsland_Grav_anom, GGM_Gi_interpolatedGippsland, 'outputs/plots/EGM2008Gippsland.png');

plotGGMGravityComparison( ...
    NSW_Long, NSW_Lat, NSW_Grav_anom, GGM_Gi_interpolatedNSW, 'outputs/plots/EGM2008NSW.png');

% Combine data sets.
ABGrav=[longitude,latitude,elev_MSL,GrvFAL100s_GEO_V,GrvFAL100s_GEO_V*0+3;...
        Gippsland_Long,Gippsland_Lat,Gippsland_H,Gippsland_Grav_anom,Gippsland_Grav_anom*0+3;...
        WA_Long, WA_Lat, WA_H, WA_Grav_anom,WA_Grav_anom*0+3;...
        Adelaide_Long,Adelaide_Lat,Adelaide_H,Adelaide_Grav_anom,Adelaide_Grav_anom*0+3;...
        NSW_Long,NSW_Lat,NSW_H,NSW_Grav_anom,NSW_Grav_anom*0+3];

% % Block mean as other datasets.
AB_Grav_Long=round(ABGrav(:,1)*60)/60;   
AB_Grav_Lat=round(ABGrav(:,2)*60)/60;
[vals, uid, idx]=unique([AB_Grav_Long,AB_Grav_Lat],'rows');
AB_Grav_BM(:,3) = accumarray(idx,ABGrav(:,3),[],@mean);
AB_Grav_BM(:,4) = accumarray(idx,ABGrav(:,4),[],@mean);
AB_Grav_BM(:,5) = accumarray(idx,ABGrav(:,5),[],@mean);
AB_Grav_BM(:,1) = AB_Grav_Long(uid);
AB_Grav_BM(:,2) = AB_Grav_Lat(uid);

AB_Grav_BM(:,6)=3;% Flag for airborne gravity
save('Data\processedData\AirborneAllMay2026.mat','AB_Grav_BM')

figure
hold on
scatter(AB_Grav_BM(:,1), AB_Grav_BM(:,2), 1, AB_Grav_BM(:,4))
colormap(jet)
cb1 = colorbar;
title(cb1,'mGal','FontSize',10)
title('Gravity anomaly')



% Helper function to replace %*f with %f at desired column
function fmtOut = replaceColumnFormat(fmtIn, colIndex, newFmt)
    tokens = strsplit(fmtIn);
    tokens{colIndex} = newFmt;
    fmtOut = strjoin(tokens, ' ');
end






