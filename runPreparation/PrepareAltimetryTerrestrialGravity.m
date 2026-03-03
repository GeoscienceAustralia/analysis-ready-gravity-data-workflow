%% Clean start
close all; clear; clc;
warning('off','all');

addpath('functions');

%% Coastline data
COAST_PARA.filename = 'Data/COASTLINE/CoastAus.mat';
Coastline = importdata(COAST_PARA.filename);

%% ------------------------------------------------------------------------
%  File definitions (TWO datasets)
% -------------------------------------------------------------------------
data(1).gravFile = fullfile('Data','GRAVITY','ALTIMETRY','sand311ausgrav.llf');
data(1).errFile  = fullfile('Data','GRAVITY','ALTIMETRY','sand311ausgrav.lle');
data(1).label    = 'Sandwell 311';

data(2).gravFile = fullfile('Data','GRAVITY','ALTIMETRY','sand331ausgrav.llf');
data(2).errFile  = fullfile('Data','GRAVITY','ALTIMETRY','sand331ausgrav.lle');
data(2).label    = 'Sandwell 331';

%% Grid dimensions
nLon = 4860;
nLat = 3180;

%% ------------------------------------------------------------------------
%  Region of interest
% -------------------------------------------------------------------------
roi.lon = [144 145];
roi.lat = [-38.5 -37.5];

GRID_PARA.buffer = 0;
GRID_PARA.MINLAT = roi.lat(1);
GRID_PARA.MAXLAT = roi.lat(2);

axisLimits.latMeanCosine = abs(cos(deg2rad(mean(roi.lat))));
axisLimits.lonMinLimit  = roi.lon(1);
axisLimits.lonMaxLimit  = roi.lon(2);
axisLimits.latMinLimit  = roi.lat(1);
axisLimits.latMaxLimit  = roi.lat(2);

%% ------------------------------------------------------------------------
%  Load, reshape, and subset both datasets
% -------------------------------------------------------------------------
for k = 1:2

    gravLLV = importdata(data(k).gravFile);
    errLLV  = importdata(data(k).errFile);

    lonGrid = reshape(gravLLV(:,1), nLon, nLat).';
    latGrid = reshape(gravLLV(:,2), nLon, nLat).';
    grav    = reshape(gravLLV(:,3), nLon, nLat).';
    err     = reshape(errLLV(:,3),  nLon, nLat).';

    lonVec = lonGrid(1,:);
    latVec = latGrid(:,1);

    iLon = lonVec >= roi.lon(1) & lonVec <= roi.lon(2);
    iLat = latVec >= roi.lat(1) & latVec <= roi.lat(2);

    data(k).lon  = lonVec(iLon);
    data(k).lat  = latVec(iLat);
    data(k).grav = grav(iLat,iLon);
    data(k).err  = err(iLat,iLon);

end

%% ------------------------------------------------------------------------
%  PLOT 1–4: gravity + uncertainty for both datasets
% -------------------------------------------------------------------------
for k = 1:2

    % ---- Uncertainty ----
    figure('Color','w');
    imagesc(data(k).lon, data(k).lat, data(k).err);
    set(gca,'YDir','normal')
    hold on
    customizeMap([data(k).label ' gravity uncertainty'], ...
                 'mGal', Coastline, axisLimits)
    saveas(gcf,[data(k).errFile,'.png'])

    % ---- Gravity ----
    figure('Color','w');
    imagesc(data(k).lon, data(k).lat, data(k).grav);
    set(gca,'YDir','normal')
    caxis([-100 100])
    hold on
    customizeMap([data(k).label ' gravity'], ...
                 'mGal', Coastline, axisLimits)
    saveas(gcf,[data(k).gravFile,'.png'])

end

%% ------------------------------------------------------------------------
%  PLOT 5–6: Differences (331 − 311)
% -------------------------------------------------------------------------
dGrav = data(2).grav - data(1).grav;
dErr  = data(2).err  - data(1).err;

% ---- Gravity difference ----
figure('Color','w');
imagesc(data(1).lon, data(1).lat, dGrav);
set(gca,'YDir','normal')
%caxis([-5 5])
hold on
customizeMap('Gravity difference (331 − 311)', ...
             'mGal', Coastline, axisLimits)
saveas(gcf,'gravity_difference_331_minus_311.png')

% ---- Uncertainty difference ----
figure('Color','w');
imagesc(data(1).lon, data(1).lat, dErr);
set(gca,'YDir','normal')
%caxis([-1 1])
hold on
customizeMap('Uncertainty difference (331 − 311)', ...
             'mGal', Coastline, axisLimits)
saveas(gcf,'uncertainty_difference_331_minus_311.png')





















