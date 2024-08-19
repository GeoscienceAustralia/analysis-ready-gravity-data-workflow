%% Gradiometry test
%
% RunGradiometryLSC is the main m-file preforming the test for gradiometry data.
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
% Main functions
% - importAndFormatData
% - computeFullTerrainEffects
% - ComputeGradiometryLSC
%
% Written by Jack McCubbine
% Last updated by Neda Darbeheshti
% Geoscience Australia, 2024-02.
%%
close all
clear all
% Turn off the irritating warnings.
warning off
%% Add the path to the function files.
addpath('functions');
plotsFolder='Outputs/plots/';
%% Grid/Tiling Parameters
% Tiling Parameters - fixed for each computation
GRID_PARA.buffer=1;% degs. The x/y extent to extract data around the tile centre.
GRID_PARA.buffer2=0.5;% degs. The x/y tile extents that are kept - where the good data are.
GRID_PARA.STEP=0.5;% The step size. This must be less than buffer2 to avoid gaps in the final grid.
GRID_PARA.filterSize=15;% filter size for spatial grid weight, this value is from experiment for tiles of one degree
GRID_PARA.filterRadius=10; % filter radius for spatial grid weight, this value is from experiment for tiles of one degree
% Grid extents - ensure these values are in GRID_PARA.STEP degree
% value increments.
% Boundary for computation
% NSW=[139.5 154 -38 -27.5];
GRID_PARA.MINLONG=139.5;%140.5; 
GRID_PARA.MAXLONG=154;%141.5; 
GRID_PARA.MINLAT=-38;%-34;   
GRID_PARA.MAXLAT=-27.5;%-29;   
%% DEM data - N.B. the dem is used to specify the grid nodes.
DEM_PARA.filename='Data\DEM\AUSDEM1m.xyz';
DEM_PARA.num_cols=4861;
DEM_PARA.num_rows=3181;
%% Gravity data
% First run ./Data/GRAVITY/XXXX/PrepareGravity_XXXXX.m
% And then /Data/GRAVITY/Combine_Gravity_Data.m
% this collates all of the gravity and position data into one matlab array.
GRAV_PARA.filename='Data\GRAVITY\Gravity.mat';%'Data/GRAVITY/AIRBORNE/WA/WA_Airborne.mat';%
GRAV_PARA.TypeB=1;% This is a Type B uncertainty value (in mGal) which is added to the uncertainty values.
GRAV_PARA.Grav_Faye_TypeB=3;
%% Gravity Gradiometry data
% Add notes here
GRAV_GRAD_PARA.filename='Data\GRAVITY_GRAD\Xcalibur_FVD_GDD.mat';
GRAV_GRAD_PARA.TypeB=1;% This is a Type B uncertainty value (in mGal/m) which is added to the uncertainty values.
GRAV_GRAD_PARA.avail=true;
%% Covariance function parameters
COV_PARA.Compute_Empircal_COV_Dec=3; % Decimation factor for empirical covariance estimates. e.g. 1 is no decimation, 2 drops 50% of the data etc. see sph_empcov for logic.
COV_PARA.Fit_Empircal_COV='auto';%'auto';% process to fit covariance N & M function values 'man' for manual to fit them on the cmd line or 'auto' , '' to just use what you supply here.
COV_PARA.FitEmpiricalCOVNSearch=[10800,1,10800,];% Start, step, stop parameter sweep values for N parameter - if auto
COV_PARA.FitEmpiricalCOVMSearch=[200,20,300];% Start, step, stop parameter sweep values for M parameter - if auto
COV_PARA.N=21600;% max Legender polynonial of cov func. 
COV_PARA.M=290;% min Legender polynonial of cov func. 
COV_PARA.width=3;% Size of precomputed cov function in degrees - must be larger the the distance between any two points on a tile. 
COV_PARA.res=30/3600; % Resolution of the covariance function
COV_PARA.COV_COMPUTED_Tilewise=true;% This recomputes the covariance function for each tile.
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
LEVELLING_PARA.Lev_eval=false;% If true, the levelling data are compared to the geoid as its computed.
LEVELLING_PARA.filename='Data\GPS_LEVELLING\Lev_NSW_NG.mat';%'Data/GPS_LEVELLING/Lev_CARS.mat';% The format of these data needs to be an array with rows [Long,Lat,h-H].
LEVELLING_PARA.Plot_Stats=false;% If true, the levelling data are compared to the geoid as its computed.
LEVELLING_PARA.Compare_To_Existing_Model=false;% If true, the levelling data are also compared to another existing geoid as its computed.
LEVELLING_PARA.Existing_Model='Data\EXISTING_GEOID_MODELS\AGQG20221120.mat';% File location of the existing model.
LEVELLING_PARA.max_diff=0.15;% Threshold for an outlier with the GNSS-levelling
%% Output
OUTPUT_PARA.Grids_name='Outputs/Grids_test/';%'Outputs/Grids_test/';%'Outputs/Grids_NSW/';%'Outputs/Grids_Bathurst/'
OUTPUT_PARA.Tiles_dir_name='Outputs/Residual_Tiles_test/';%'Outputs/Residual_Tiles_test/'; %Outputs/Residual_Tiles_NSW/';%'Outputs/Residual_Tiles_Bathurst/'
PLOT_GRIDS=true;% A gridded solution is plotted and output as well as the tiles.
%% Keep the computer awake
keepawake=true;% Setting this to true wiggles the mouse every so often so the compute doesnt go to sleep.

disp('1/3 ImportAndFormatData is running ')
disp(('...................................................................'))

[Gravo,Grav_grad,DEM_data,ZDEM_griddedInterpolant,LongDEM,LatDEM,...
 GGM_Gravity_griddedInterpolant,GGM_Zeta_griddedInterpolant,Lev,...
 REFERENCE_Zeta_griddedInterpolant,GRID_REF,Coastline]=importAndFormatData...
 (GRID_PARA,DEM_PARA,GRAV_PARA,COAST_PARA,LEVELLING_PARA,GGM_PARA,GRAV_GRAD_PARA);

Nanfilter = createNanFilter(Grav_grad(:,4),470,606);
Grav_gradBlockMedian = Grav_grad(~isnan(Grav_grad(:,4).*Nanfilter), :);
% plot scatter input data

plotCustomScatter(DEM_data(:,1),DEM_data(:,2),DEM_data(:,3),GRID_PARA,Topo_PARA.Rad,'DEM','m',plotsFolder)
% hold on
% xline(Gravo(13676,1),'--',{'profile'}) 
plotCustomScatter(Gravo(:,1),Gravo(:,2),Gravo(:,3),GRID_PARA,Topo_PARA.Rad,'GravityTopographyHeight','m',plotsFolder)
% hold on
% xline(Gravo(13676,1),'--',{'profile'}) 
plotCustomScatter(Gravo(:,1),Gravo(:,2),Gravo(:,4),GRID_PARA,Topo_PARA.Rad,'Gravity','mGal',plotsFolder)
% hold on
% xline(Gravo(13676,1),'--',{'profile'}) 
plotCustomScatter(Grav_gradBlockMedian(:,1),Grav_gradBlockMedian(:,2),Grav_gradBlockMedian(:,3),GRID_PARA,Topo_PARA.Rad,'GravityGradientFlightAltitude','m',plotsFolder)
% hold on
% xline(Grav_grad(1,1),'--',{'profile'}) 
plotCustomScatter(Grav_gradBlockMedian(:,1),Grav_gradBlockMedian(:,2),Grav_gradBlockMedian(:,4),GRID_PARA,Topo_PARA.Rad,'GravityGradient','mGal/m',plotsFolder)
% hold on
% xline(Grav_grad(1,1),'--',{'profile'})
% plot two profiles input data

plotProfiles(Grav_gradBlockMedian(1:303,2),Grav_gradBlockMedian(1:303,3),Grav_gradBlockMedian(1:303,4),Grav_gradBlockMedian(1,1),'Latitude','Elevation [m]','Gradient [mGal/m]','GradiometryElevation',plotsFolder)

plotProfiles(Grav_gradBlockMedian(1:303,2),ZDEM_griddedInterpolant(Grav_gradBlockMedian(1:303,1),Grav_gradBlockMedian(1:303,2)),Grav_gradBlockMedian(1:303,4),Grav_gradBlockMedian(1,1),'Latitude','DEM [m]','Gradient [mGal/m]','GradiometryDEM',plotsFolder)


% plotProfiles(Gravo(13676:13796,2),Gravo(13676:13796,3),Gravo(13676:13796,4),Gravo(13676,1),'Latitude','Height [m]','Gravity [mGal]','GravityHeight',plotsFolder)
% 
% plotProfiles(Gravo(13676:13796,2),ZDEM_griddedInterpolant(Gravo(13676:13796,1),Gravo(13676:13796,2)),Gravo(13676:13796,4),Gravo(13676,1),'Latitude','DEM [m]','Gravity [mGal]','GravityDEM',plotsFolder)


%plotProfiles(Grav_grad(1:303,2),GGM_Gravity_griddedInterpolant(Grav_grad(1:303,1),Grav_grad(1:303,2)),Grav_grad(1:303,4),'GravityGGM [mGal]','Gradient [mGal/m]','Latitude','Gradiometry',plotsFolder)

%plotProfiles(Grav_grad(1:303,2),GGM_Zeta_griddedInterpolant(Grav_grad(1:303,1),Grav_grad(1:303,2)),Grav_grad(1:303,4),'ZetaGGM [m]','Gradient [mGal/m]','Latitude','Gradiometry',plotsFolder)

% compute 2nd order Free-Air Correction for WGS-84 Ellipsoid
%secondOrderFreeAirCorrection=computeSecondOrderFreeAirCorrection(Grav_grad(:,2),Grav_grad(:,3));

%plotCustomScatter(Grav_grad(:,1),Grav_grad(:,2),secondOrderFreeAirCorrection,GRID_PARA,'secondOrderFreeAirCorrection','mGal',plotsFolder)


% disp('2/3 computeFullTerrainEffects is running')
% disp(('...................................................................'))
% 
% [fullTopoCorrectedGravityPoint,longwaveTopo_griddedInterpolant,fullTopo_griddedInterpolant,TCgradientDEMpoint]=computeFullTerrainEffects(GRID_PARA, ...
%     Topo_PARA,Gravo,GGM_Gravity_griddedInterpolant,DEM_data,ZDEM_griddedInterpolant, ...
%     LongDEM,LatDEM,plotsFolder);


disp('3/3 ComputeGradiometryLSC is running')
disp(('...................................................................'))

computeGradiometryLSC(GRID_PARA,COV_PARA,...
    PLOT_GRIDS,Gravo,Grav_grad,Grav_gradBlockMedian, ...
    GGM_Gravity_griddedInterpolant,ZDEM_griddedInterpolant,Nanfilter,plotsFolder)

% disp('3/4 ..........................computeGradiometryLSC is running')
% computeGradiometryLSC(GRID_PARA,COV_PARA,DEM_PARA,GRAV_PARA,GRAV_GRAD_PARA,OUTPUT_PARA,GRID_REF,TE.fullTopoCorrectedGravityPoint,TE.fullTopoCorrectedGravityGradient, ...
%     GGM_Gravity_griddedInterpolant,ZDEM_griddedInterpolant,TE.fullTopo_griddedInterpolant, ...
%     TE.longwaveTopo_griddedInterpolant,Topo_PARA.Density)


