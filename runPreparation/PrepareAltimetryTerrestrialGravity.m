%% Clean start
close all; clear; clc;
warning('off','all');  % (optional) avoid blanket warning off if you prefer

addpath('functions');
% Coastline data
COAST_PARA.filename='Data/COASTLINE/CoastAus.mat';
%% ------------------------------------------------------------------------
%  Load Sandwell & Smith altimetry free-air gravity + uncertainty
%  Files appear to be [lon lat value] columns
% -------------------------------------------------------------------------
%gravFile = fullfile('Data','GRAVITY','ALTIMETRY','sand311ausgrav.llf');
%errFile  = fullfile('Data','GRAVITY','ALTIMETRY','sand311ausgrav.lle');

 gravFile = fullfile('Data','GRAVITY','ALTIMETRY','sand331ausgrav.llf');
 errFile  = fullfile('Data','GRAVITY','ALTIMETRY','sand331ausgrav.lle');

gravLLV = importdata(gravFile);   % [lon lat freeAir]
errLLV  = importdata(errFile);    % [lon lat freeAirError]

% Grid dimensions used in your reshape
nLon = 4860;
nLat = 3180;

% Reshape into 2-D grids (lat x lon) — matches your transpose pattern
lonGrid      = reshape(gravLLV(:,1), nLon, nLat).';
latGrid      = reshape(gravLLV(:,2), nLon, nLat).';
freeAir_mGal = reshape(gravLLV(:,3), nLon, nLat).';

freeAirErr_mGal = reshape(errLLV(:,3), nLon, nLat).';

%% ------------------------------------------------------------------------
%  Choose region of interest: South Australia (adjust as needed)
% -------------------------------------------------------------------------
% Broad South Australia-ish bounding box (degrees)
roi.lon = [144 145];
roi.lat = [-38.5 -37.5];

% If you want a tighter "Greater Adelaide" view, use instead:
% roi.lon = [137.5 139.5];
% roi.lat = [-35.8 -34.2];

% Build mask / index ranges using the grids
lonVec = lonGrid(1,:);   % longitudes along columns
latVec = latGrid(:,1);   % latitudes along rows

iLon = lonVec >= roi.lon(1) & lonVec <= roi.lon(2);
iLat = latVec >= roi.lat(1) & latVec <= roi.lat(2);

lonSub = lonVec(iLon);
latSub = latVec(iLat);
errSub = freeAirErr_mGal(iLat, iLon);
gravSub = freeAir_mGal(iLat, iLon);

% common variables for plotting
Coastline=importdata(COAST_PARA.filename);
GRID_PARA.buffer=0;
GRID_PARA.MINLAT=roi.lat(1);
GRID_PARA.MAXLAT=roi.lat(2);

axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
axisLimits.lonMinLimit=roi.lon(1)-GRID_PARA.buffer;
axisLimits.lonMaxLimit=roi.lon(2)+GRID_PARA.buffer;
axisLimits.latMinLimit=roi.lat(1)-GRID_PARA.buffer;
axisLimits.latMaxLimit=roi.lat(2)+GRID_PARA.buffer;


figure('Color','w');
imagesc(lonSub, latSub, errSub);
set(gca,'YDir','normal')
hold on
customizeMap('gravity uncertainty','mGal',Coastline,axisLimits)
saveas(gcf,[errFile,'.png'])

figure('Color','w');
imagesc(lonSub, latSub, gravSub);
set(gca,'YDir','normal')
caxis([-100 100]);
hold on
customizeMap('gravity','mGal',Coastline,axisLimits)
saveas(gcf,[gravFile,'.png'])






















