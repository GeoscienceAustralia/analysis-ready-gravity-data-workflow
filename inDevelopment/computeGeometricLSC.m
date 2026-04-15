%% Gravimetric quasigeoid estimation
%
% RunGeoidLSC is the main m-file preforming the gravimetric quasiGeoid estimation.
%
% Main parameters to be defined:
%                      - Grid/Tiling Parameters
%                      - DEM (Digital Elevation Model) data
%                      - Gravity data
%                      - Gravity Gradiometry data
%                      - Covariance function parameters
%                      - Topo condensation parameters
%                      - GGM (Global Gravity Model) reference signal
%                      - Coastline data
%                      - Levelling data comparisons
%                      - Output directories
%
% Main functions:
%             - ImportAndFormatData: import all different data sets needed for LSC.
%             - computeTerrainEffect: calculates the Residual Terrain Model.
%             - ComputeLSC: runs the LSC (Least Squares Collocation) in blocks.
%             - MosaicTiles: mosaic the tiles with weights and save final result in tif 
%
% Written by Jack McCubbine
% Last updated by Neda Darbeheshti
% Geoscience Australia, 2024-03.
%%
close all
clear 
% Turn off the irritating warnings.
warning off
%% Add the path to the function files.
addpath('functions');
%% Grid/Tiling Parameters
% Tiling Parameters - fixed for each computation
GRID_PARA.buffer=1;% degs. The x/y extent to extract data around the tile centre.
GRID_PARA.buffer2=0.5;% degs. The x/y tile extents that are kept - where the good data are.
GRID_PARA.STEP=0.5;% The step size. This must be less than buffer2 to avoid gaps in the final grid.
GRID_PARA.filterSize=15;% filter size for spatial grid weight, this value is from experiment for tiles of one degree
GRID_PARA.filterRadius=10; % filter radius for spatial grid weight, this value is from experiment for tiles of one degree
% Grid extents - ensure these values are in GRID_PARA.STEP degree value increments.
% Boundary for computation
%Aus=[114 154 -44 -10];
GRID_PARA.MINLONG=114;
GRID_PARA.MAXLONG=154;
GRID_PARA.MINLAT=-44;
GRID_PARA.MAXLAT=-10;
%% DEM data - N.B. the dem is used to specify the grid nodes.
DEM_PARA.filename='Data\DEM\AUSDEM1min.xyz';
DEM_PARA.num_cols=4861;
DEM_PARA.num_rows=3181;
%% Gravity data
% First run ./Data/GRAVITY/XXXX/PrepareGravity_XXXXX.m
% And then /Data/GRAVITY/Combine_Gravity_Data.m
% this collates all of the gravity and position data into one matlab array.
GRAV_PARA.filename='Data\processedData\GravityAllTerrestrialAirborneJuly14.mat';
GRAV_PARA.filename1=[];%'Data\GRAVITY\Xcalibur_Gravity.mat';
GRAV_PARA.TypeB=1;% This is a Type B uncertainty value (in mGal) which is added to the uncertainty values.
GRAV_PARA.Grav_Faye_TypeB=3;
%% Gravity Gradiometry data
% Add notes here
GRAV_GRAD_PARA.filename='Data/processedData/OtwayXcalibur.mat';
GRAV_GRAD_PARA.TypeB=10^(-4);% This is a Type B uncertainty value (in mGal/m) which is added to the uncertainty values.
GRAV_GRAD_PARA.avail=true;
%% Covariance function parameters
COV_PARA.Compute_Empircal_COV_Dec=3; % Decimation factor for empirical covariance estimates. e.g. 1 is no decimation, 2 drops 50% of the data etc. see sph_empcov for logic.
COV_PARA.Fit_Empircal_COV='auto';%'auto';% process to fit covariance N & M function values 'man' for manual to fit them on the cmd line or 'auto' , '' to just use what you supply here.
COV_PARA.FitEmpiricalCOVNSearch=[21600,1,21600]; %21600% Start, step, stop parameter sweep values for N parameter - if auto
COV_PARA.FitEmpiricalCOVMSearch=[200,20,300];% Start, step, stop parameter sweep values for M parameter - if auto
COV_PARA.N=10800;% max Legender polynonial of cov func. 
COV_PARA.M=200;% min Legender polynonial of cov func. 
COV_PARA.width=3;% Size of precomputed cov function in degrees - must be larger the the distance between any two points on a tile. 
COV_PARA.res=30/3600; % Resolution of the covariance function
COV_PARA.COV_COMPUTED_Tilewise=true; %false% This recomputes the covariance function for each tile.
COV_PARA.Airbornedataonly=false;%Only use airborne data in establishing Covariance parameters - good to use if we are using EGM2008 as the references as terrestrial data are not independent.
COV_PARA.COVPlot=true;% true plots progress, false turns this off.
%% Topo condensation parameters
Topo_PARA.Corr=true;% MAKE SURE YOU TURN THIS ON!!!
Topo_PARA.TopoPlot=true;% true plots progress, false turns this off.
Topo_PARA.Density=2.67;% Assumed density in g/cm^3.
Topo_PARA.Depth=0;% Condensation layer depth. 0 is on the geoid
Topo_PARA.Rad=1;% Radius out to which to compute the effects in degress.
Topo_PARA.RTM=[50,10,300];%[1000,10,2160]for egm%[0,10,300]% Range of SHM degree filter parameters (min, step, max) explored when running RTM calculations.
%about 1080 for EGM,we need to fix the filter by factor of 2 
%% GGM reference signal
% First run e.g. ./Data/GGM/RunIsGrafLab_Topo_Surf_EGM2008.m
GGM_PARA.filename='Data\GGM\GOCE_For_Gridded_Int.mat';%'Data/GGM/EGM2008_For_Gridded_Int.mat';%'Data\GGM\GOCE_For_Gridded_Int.mat';%'Data/GGM/GOCE_N200_For_Gridded_Int.mat';%'Data/GGM/Tongji_For_Gridded_Int.mat';%'Data/GGM/XGM2019_For_Gridded_Int.mat';
%% Coastline data
COAST_PARA.filename='Data\COASTLINE\CoastAus.mat';
%% Levelling data comparisons
LEVELLING_PARA.Lev_eval=true;% If true, the levelling data are compared to the geoid as its computed.
LEVELLING_PARA.filename='Data/GPS_LEVELLING/AHDzeta7319.mat';%'Data/GPS_LEVELLING/Lev_CARS.mat';% The format of these data needs to be an array with rows [Long,Lat,h-H].
LEVELLING_PARA.Plot_Stats=true;% If true, the levelling data are compared to the geoid as its computed.
LEVELLING_PARA.Compare_To_Existing_Model=true;% If true, the levelling data are also compared to another existing geoid as its computed.
LEVELLING_PARA.Existing_Model='Data\EXISTING_GEOID_MODELS\AGQG20221120.mat';% File location of the existing model.
LEVELLING_PARA.max_diff=0.15;% Threshold for an outlier with the GNSS-levelling
%% Output
outputName='AustraliaSparse';
plotName='ahd';
OUTPUT_PARA.Grids_name=['outputs/Grids',outputName,'/'];
OUTPUT_PARA.PLOT_GRIDS=true;% A gridded solution is plotted and output as well as the tiles.
OUTPUT_PARA.plotsFolder=['outputs/Grids',outputName,'/',plotName];
% Keep the computer awake
keepawake=true;% Setting this to true wiggles the mouse every so often so the compute doesnt go to sleep.

disp('1/4 ..........................importAndFormatData is running ')
[Gravo,gravGradFiltered,DEM_data,ZDEM_griddedInterpolant,LongDEM,LatDEM,...
 GGM_Gravity_griddedInterpolant,GGM_Zeta_griddedInterpolant,Lev,...
 REFERENCE_Zeta_griddedInterpolant,GRID_REF,Coastline]=importAndFormatData...
 (GRID_PARA,DEM_PARA,GRAV_PARA,Topo_PARA,COAST_PARA,LEVELLING_PARA,GGM_PARA,GRAV_GRAD_PARA);

% read final matfiles

dateCreated ='23-Mar-2026';

load([OUTPUT_PARA.Grids_name,'Grid_res_geoid_w',dateCreated,'.mat'])

ZDeg=mean(mean(REFERENCE_Zeta_griddedInterpolant(LongDEM,LatDEM)-GGM_Zeta_griddedInterpolant(LongDEM,-LatDEM,LatDEM*0)));

resAGQG=REFERENCE_Zeta_griddedInterpolant(LongDEM,LatDEM)-GGM_Zeta_griddedInterpolant(LongDEM,-LatDEM,LatDEM*0);

Geoid_temp=double(Grid_res_geoid_w+GGM_Zeta_griddedInterpolant(LongDEM,-LatDEM,LatDEM*0));
       
geoidLSCgriddedInterpolant=griddedInterpolant(LongDEM(end:-1:1,:)',LatDEM(end:-1:1,:)',Geoid_temp(end:-1:1,:)');
    
geomGravDiff=Lev(:,3)-geoidLSCgriddedInterpolant(Lev(:,1),Lev(:,2));  

geomGravDiff2022=Lev(:,3)-REFERENCE_Zeta_griddedInterpolant(Lev(:,1),Lev(:,2)); 

% plot
plotGPSlevelling(Coastline,GRID_PARA,Lev,geomGravDiff,geomGravDiff2022,OUTPUT_PARA.plotsFolder)

% Remove a tiled plane so the signal is zero mean for the LSC

% Construct the matrix for linear trend removal
trendMatrix = [Lev(:,1) - mean(Lev(:,1)), Lev(:,2) - mean(Lev(:,2)), ones(size(Lev(:,2)))];

% Calculate the coefficients of the best-fit plane
trendCoefficients = trendMatrix \ geomGravDiff;
trendCoefficients2022 = trendMatrix \ geomGravDiff2022;

% Remove the planar trend to obtain zero-mean data
geomGravGeoidDiffDetrended = geomGravDiff - trendMatrix * trendCoefficients;
geomGravGeoidDiff2022Detrended = geomGravDiff2022 - trendMatrix * trendCoefficients2022;

% plot
plotGPSlevelling(Coastline,GRID_PARA,Lev,geomGravGeoidDiffDetrended,geomGravGeoidDiff2022Detrended,OUTPUT_PARA.plotsFolder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
disp('computing covariance functions')

covarianceInfo=computeSphericalEmpiricalCovariance(Lev(:,1),Lev(:,2),geomGravGeoidDiffDetrended,1);

[sigma2,bestFitCoeff,fittedCovariance]=fitGaussianCovariance(covarianceInfo(:,1),covarianceInfo(:,2));

% Plot covariance function
plotSphericalCovarianceFunction(covarianceInfo(:,1), covarianceInfo(:,2), fittedCovariance,'m^2','global Gaussian Covariance', OUTPUT_PARA.plotsFolder)

% Convert degrees to radians
longitudeLevRadian = deg2rad (Lev(:,1));
latitudeLevRadian = deg2rad (Lev(:,2));

% initialize covariance matrix
ACOVtt = zeros(length(Lev(:,1)),length(Lev(:,1)));
for lonCounter=1:length(Lev(:,1))
haversineDistance=haversine(latitudeLevRadian(lonCounter), longitudeLevRadian(lonCounter),latitudeLevRadian(:), longitudeLevRadian(:));
ACOVtt(lonCounter,:)=sigma2*exp(-(haversineDistance.^2)/(2*bestFitCoeff.^2));
end

% Plot covariance function
plotSphericalCovarianceFunction(haversineDistance, ACOVtt(lonCounter, :), 0*rad2deg(haversineDistance),'m^2','ACOVtt', OUTPUT_PARA.plotsFolder)

% LSC matrix multiplication 
% inverse of auto covariance matrix
inverseCovarianceMatrix=(ACOVtt+0.000025*eye(size(ACOVtt)))\eye(size(ACOVtt));
 
temporaryVector=inverseCovarianceMatrix*(geomGravGeoidDiffDetrended);
%%%%%%%%%%%%% this block trimmes and cut the DEM
disp('DEM')
DEM3D=importdata(DEM_PARA.filename);
disp('Extracting DEM subset') 
%make sure DEM is bigger than gravity
Topo_buffer=Topo_PARA.Rad+GRID_PARA.buffer; 
CoordsMM_topo=[GRID_PARA.MINLONG-Topo_buffer,GRID_PARA.MINLAT-Topo_buffer;...
          GRID_PARA.MINLONG-Topo_buffer,GRID_PARA.MAXLAT+Topo_buffer;...
          GRID_PARA.MAXLONG+Topo_buffer,GRID_PARA.MAXLAT+Topo_buffer;...
          GRID_PARA.MAXLONG+Topo_buffer,GRID_PARA.MINLAT-Topo_buffer;...
          GRID_PARA.MINLONG-Topo_buffer,GRID_PARA.MINLAT-Topo_buffer];

DEMin=inpolygon(DEM3D(:,1),DEM3D(:,2),CoordsMM_topo(:,1),CoordsMM_topo(:,2));
DEM3D(DEMin==0,:)=[];
%Computing grid dimensions for one-minute spatial resolution
DEM_PARA.num_cols=(max(DEM3D(:,1))-min(DEM3D(:,1)))*60+1;
DEM_PARA.num_rows=(max(DEM3D(:,2))-min(DEM3D(:,2)))*60+1;
%Set the computational grid nodes
LongDEM=reshape(DEM3D(:,1),DEM_PARA.num_cols,DEM_PARA.num_rows)';
LatDEM=reshape(DEM3D(:,2),DEM_PARA.num_cols,DEM_PARA.num_rows)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Extract the datasets over the tile region')
Topo_buffer=Topo_PARA.Rad+GRID_PARA.buffer; 
CoordsMM_topo=[GRID_PARA.MINLONG-Topo_buffer,GRID_PARA.MINLAT-Topo_buffer;...
          GRID_PARA.MINLONG-Topo_buffer,GRID_PARA.MAXLAT+Topo_buffer;...
          GRID_PARA.MAXLONG+Topo_buffer,GRID_PARA.MAXLAT+Topo_buffer;...
          GRID_PARA.MAXLONG+Topo_buffer,GRID_PARA.MINLAT-Topo_buffer;...
          GRID_PARA.MINLONG-Topo_buffer,GRID_PARA.MINLAT-Topo_buffer];


INOUT=inpolygon(GRID_REF(:,1),GRID_REF(:,2),CoordsMM_topo(:,1),CoordsMM_topo(:,2)); % Whole tile zone mask - n.b. data on edges of tile are unreliable
GRID_REF_dat=GRID_REF(INOUT==1,:);
GRID_REF_dat(:,1)=round(GRID_REF_dat(:,1)*60)/60;
GRID_REF_dat(:,2)=round(GRID_REF_dat(:,2)*60)/60;
%%%%%%%%%%%%%%%%%%%%%%%% new way to make the covariance matrix
constants                                       % load constants
phi=deg2rad(mean(GRID_REF(:,2)));
RadiusBjerhammar= EarthMajorAxis*EarthMinorAxis/sqrt((EarthMajorAxis*sin(phi)).^2+(EarthMinorAxis*cos(phi)).^2)*10^3;% Pajama sphere radius.

CCov_tt_int_fun_RTM=precomputeCovarianceFunction('cov_tt',RadiusBjerhammar,COV_PARA.width,COV_PARA.res,sigma2,bestFitCoeff,COV_PARA.N,COV_PARA.M);

% Auto-covariance of potential at DEM points
ACOVttRTM_dem = interpolateCovarianceFunction(...
Lev(:,1), Lev(:,2), ...
RadiusBjerhammar + ZDEM_griddedInterpolant(Lev(:,1), Lev(:,2)), ...
Lev(:,1), Lev(:,2), ...
RadiusBjerhammar + ZDEM_griddedInterpolant(Lev(:,1), Lev(:,2)), CCov_tt_int_fun_RTM,OUTPUT_PARA,'ACOVttGPSlevelling m^4/s^4',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%doing the multiplication one row of latitude at a time.
%Convert degrees to radians
longitudeLongDEMRadian = deg2rad (LongDEM);
latitudeLatDEMRadian = deg2rad (LatDEM);

ACOV_tt_dem = zeros(size(LongDEM,2),length(Lev(:,1)));
   
LSC_sol=LongDEM*0;
LSC_solrt=LSC_sol;

%for latCounter=1:length(LongDEM(:,1))
 for latCounter=1:1

    ACOV_tt_dem=[];

    for lonCounter=1:length(Lev(:,1))

    haversineDistance=haversine(latitudeLevRadian(lonCounter), longitudeLevRadian(lonCounter),latitudeLatDEMRadian(latCounter,:), longitudeLongDEMRadian(latCounter,:));
    ACOV_tt_dem(lonCounter,:)=sigma2*exp(-(haversineDistance.^2)/(2*bestFitCoeff.^2));
    
    end

    ACOV_tt_dem=ACOV_tt_dem';     
    LSC_sol(latCounter,:)=ACOV_tt_dem*temporaryVector;
    disp(latCounter)
    % Add code to restore the tilt
    LSC_solrt(:)=LSC_sol(:)+[LongDEM(:)-mean(Lev(:,1)),LatDEM(:)-mean(Lev(:,2)),ones(size(LongDEM(:)))]*trendCoefficients;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
longitudeLongDEMRadian = deg2rad (LongDEM);
latitudeLatDEMRadian = deg2rad (LatDEM);

ACOV_tt_dem = zeros(size(LongDEM,2),length(Lev(:,1)));
   
LSC_sol=LongDEM*0;
LSC_solrt=LSC_sol;


nLat = size(LongDEM, 1);

for latCounter = 1:1
    % Vectorized haversine over all Lev points at once (if haversine supports it)
    D = haversine( ...
        latitudeLevRadian(:), longitudeLevRadian(:), ...
        latitudeLatDEMRadian(latCounter,:), longitudeLongDEMRadian(latCounter,:) );

    % Build covariance in one go
    ACOV_tt_dem = sigma2 * exp(-(D.^2) / (2*bestFitCoeff.^2));

    % Ensure correct orientation (depends on haversine output shape)
    LSC_sol(latCounter,:) = (ACOV_tt_dem.') * temporaryVector;

    disp(latCounter)
end


parfor latCounter = 1:nLat

    % Vectorised haversine over all Lev points
    D = haversine( ...
        latitudeLevRadian(:), longitudeLevRadian(:), ...
        latitudeLatDEMRadian(latCounter,:), longitudeLongDEMRadian(latCounter,:) );

    % Build covariance
    ACOV_tt_dem = sigma2 * exp(-(D.^2) / (2*bestFitCoeff.^2));

    % Store result (sliced output – parfor safe)
    LSC_sol(latCounter,:) = (ACOV_tt_dem.') * temporaryVector;

end

LSC_solrt = LSC_sol(:) + ...
    [LongDEM(:) - mean(Lev(:,1)), ...
     LatDEM(:)  - mean(Lev(:,2)), ...
     ones(size(LongDEM(:)))] * trendCoefficients;



save([OUTPUT_PARA.Grids_name,'geometricgeoidgg',date,'.mat'],'LongDEM','LatDEM','LSC_solrt','LSC_sol')

fprintf('%f size    gridgeomGravGeoidDiffDetrended\n',size          (LSC_sol));
fprintf('%f min     gridgeomGravGeoidDiffDetrended\n',min(min       (LSC_sol)));
fprintf('%f max     gridgeomGravGeoidDiffDetrended\n',max(max       (LSC_sol)));
fprintf('%f mean    gridgeomGravGeoidDiffDetrended\n',mean(mean     (LSC_sol)));
fprintf('%f median  gridgeomGravGeoidDiffDetrended\n',median(median (LSC_sol)));
fprintf('%f std     gridgeomGravGeoidDiffDetrended\n',std(std       (LSC_sol)));

% common variables for plotting
axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;

figure('Name','Grid','NumberTitle','off'); 
clf
hold on
imagesc(LongDEM(1,:),LatDEM(:,1),LSC_sol)
customizeMap('geomGravGeoidDiffDetrended','m',Coastline,axisLimits)
caxis([-.1 .1])
saveas(gcf,[OUTPUT_PARA.plotsFolder,'Grid','geomGravGeoidDiffDetrended','.png'])

figure('Name','Grid','NumberTitle','off'); 
clf
hold on
imagesc(LongDEM(1,:),LatDEM(:,1),LSC_solrt)
customizeMap('geomGravGeoidDiff','m',Coastline,axisLimits)
saveas(gcf,[OUTPUT_PARA.plotsFolder,'Grid','geomGravGeoidDiff','.png'])




 



