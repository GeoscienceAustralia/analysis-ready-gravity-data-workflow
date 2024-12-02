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
GRID_PARA.buffer=1;% 1 degs. The x/y extent to extract data around the tile centre. .75
GRID_PARA.buffer2=0.5;% degs. The x/y tile extents that are kept - where the good data are.
GRID_PARA.STEP=0.5;% The step size. This must be less than buffer2 to avoid gaps in the final grid.
GRID_PARA.filterSize=15;% filter size for spatial grid weight, this value is from experiment for tiles of one degree
GRID_PARA.filterRadius=10; % filter radius for spatial grid weight, this value is from experiment for tiles of one degree
% Grid extents - ensure these values are in GRID_PARA.STEP degree value increments.
% Boundary for computation
% VicNSW=[140 154 -38 -27.5];
% NENSW=[153 154 -29 -28];
GRID_PARA.MINLONG=115;%153%115
GRID_PARA.MAXLONG=116;%154%117
GRID_PARA.MINLAT=-33;%-29%-33
GRID_PARA.MAXLAT=-32;%-28%-31
%% DEM data - N.B. the dem is used to specify the grid nodes.
DEM_PARA.filename='Data/DEM/AUSDEM1min.xyz';
DEM_PARA.num_cols=4861;
DEM_PARA.num_rows=3181;
%% Gravity data
% First run ./Data/GRAVITY/XXXX/PrepareGravity_XXXXX.m
% And then /Data/GRAVITY/Combine_Gravity_Data.m
% this collates all of the gravity and position data into one matlab array.
GRAV_PARA.filename = 'Data/processedData/GravityAllVicNSW.mat';
GRAV_PARA.filename1 = [];%'Data/GRAVITY/Xcalibur_Gravity.mat';% gravity from gradiometry
GRAV_PARA.TypeB = 10;%1;% This is a Type B uncertainty value (in mGal) which is added to the uncertainty values.
GRAV_PARA.Grav_Faye_TypeB = 30;%3;
GRAV_PARA.inputGravity_weighting = false; 
%% Gravity Gradiometry data
% Add notes here
GRAV_GRAD_PARA.filename='Data/GRAVITY_GRAD/Xcalibur_FVD_GDD.mat';
GRAV_GRAD_PARA.TypeB=10^(-5);% This is a Type B uncertainty value (in mGal/m) which is added to the uncertainty values.
GRAV_GRAD_PARA.avail=true;
%% Covariance function parameters
COV_PARA.Compute_Empircal_COV_Dec=3; % Decimation factor for empirical covariance estimates. e.g. 1 is no decimation, 2 drops 50% of the data etc. see sph_empcov for logic.
COV_PARA.Fit_Empircal_COV='auto';%'auto';% process to fit covariance N & M function values 'man' for manual to fit them on the cmd line or 'auto' , '' to just use what you supply here.
COV_PARA.FitEmpiricalCOVNSearch=[21600,1,21600];%10800 is not working:speckelly pattern% Start, step, stop parameter sweep values for N parameter - if auto
COV_PARA.FitEmpiricalCOVMSearch=[200,20,300];% Start, step, stop parameter sweep values for M parameter - if auto
COV_PARA.N=10800;% max Legender polynonial of cov func. 
COV_PARA.M=200;% min Legender polynonial of cov func. 
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
GGM_PARA.filename='Data/GGM/GOCE_For_Gridded_Int.mat';%'Data/GGM/EGM2008_For_Gridded_Int.mat';%'Data/GGM/GOCE_For_Gridded_Int.mat';%'Data/GGM/GOCE_N200_For_Gridded_Int.mat';%'Data/GGM/Tongji_For_Gridded_Int.mat';%'Data/GGM/XGM2019_For_Gridded_Int.mat';
%% Coastline data
COAST_PARA.filename='Data/COASTLINE/CoastAus.mat';
%% Levelling data comparisons
LEVELLING_PARA.Lev_eval=true;% If true, the levelling data are compared to the geoid as its computed.
LEVELLING_PARA.filename='Data/GPS_LEVELLING/Lev_NSW_NG.mat';%'Data/GPS_LEVELLING/Lev_CARS.mat';% The format of these data needs to be an array with rows [Long,Lat,h-H].
LEVELLING_PARA.Plot_Stats=false;% If true, the levelling data are compared to the geoid as its computed.
LEVELLING_PARA.Compare_To_Existing_Model=true;% If true, the levelling data are also compared to another existing geoid as its computed.
LEVELLING_PARA.Existing_Model='Data/EXISTING_GEOID_MODELS/AGQG20221120.mat';% File location of the existing model.
LEVELLING_PARA.max_diff=0.15;% Threshold for an outlier with the GNSS-levelling
%% Output
OUTPUT_PARA.Grids_name='outputs/GridsPerthSim165altTypeB1030/';
OUTPUT_PARA.Tiles_dir_name='outputs/ResidualTilesPerthSim165altTypeB1030/';
OUTPUT_PARA.PLOT_GRIDS=true;% A gridded solution is plotted and output as well as the tiles.
OUTPUT_PARA.plotsFolder=['outputs/plots/',date,'PerthSim165altTypeB1030'];
% Keep the computer awake
keepawake=true;% Setting this to true wiggles the mouse every so often so the compute doesnt go to sleep.

diary([OUTPUT_PARA.Grids_name,date,'PerthSim165altTypeB1030.txt']);% start recording  

disp('1/4 ..........................importAndFormatData is running ')
[Gravo,gravGradFiltered,DEM_data,ZDEM_griddedInterpolant,LongDEM,LatDEM,...
 GGM_Gravity_griddedInterpolant,GGM_Zeta_griddedInterpolant,Lev,...
 REFERENCE_Zeta_griddedInterpolant,GRID_REF,Coastline,DEM_PARA]=importAndFormatData...
 (GRID_PARA,DEM_PARA,GRAV_PARA,Topo_PARA,COAST_PARA,LEVELLING_PARA,GGM_PARA,GRAV_GRAD_PARA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% create simulated airborne gravity data
% Assuming your data is in a matrix format with columns: longitude, latitude, height, gravity
data = Gravo; % Replace with your actual data file

% Extract columns
lon = data(:, 1);
lat = data(:, 2);
height = data(:, 3);
gravity = data(:, 4);

% Define the grid for interpolation (1-minute arc)
lon_grid = min(lon):1/60:max(lon);
lat_grid = min(lat):1/60:max(lat);
[LonGrid, LatGrid] = meshgrid(lon_grid, lat_grid);

% Interpolate gravity data onto the grid
GravityGrid = griddata(lon, lat, gravity, LonGrid, LatGrid, 'natural');

% Interpolate height data onto the grid
HeightGrid = griddata(lon, lat, height, LonGrid, LatGrid, 'natural');

% Flatten the grids to match the original data format
LonGrid_flat = LonGrid(:);
LatGrid_flat = LatGrid(:);
GravityGrid_flat = GravityGrid(:);
HeightGrid_flat = HeightGrid(:);

% Create two columns with the number 3
column_airborne_uncertainty = ones(size(LonGrid_flat)) * 3;
column_airborne_flag = ones(size(LonGrid_flat)) * 3;

% Combine the gridded data into a single matrix with the additional columns
Gravo = [LonGrid_flat, LatGrid_flat, HeightGrid_flat+165, GravityGrid_flat, column_airborne_uncertainty, column_airborne_flag];

if OUTPUT_PARA.PLOT_GRIDS
     plotInputData(Gravo,Coastline,GRID_PARA,OUTPUT_PARA)
end 

if GRAV_PARA.inputGravity_weighting 
     Gravo = weightInputGravity(Gravo,Coastline,GRID_PARA,OUTPUT_PARA);
end

if exist([OUTPUT_PARA.Grids_name,'terrainEffects.mat'], 'file')
    load([OUTPUT_PARA.Grids_name,'terrainEffects.mat']);
    disp('3/4 ..........................computeGravimetryGradiometryLSC is running')
    computeGravimetryGradiometryLSC(GRID_PARA,COV_PARA,DEM_PARA,GRAV_PARA,GRAV_GRAD_PARA,OUTPUT_PARA,GRID_REF,fullTopoCorrectedGravityPoint,fullTopoCorrectedGravityGradient, ...
        GGM_Gravity_griddedInterpolant,ZDEM_griddedInterpolant,fullTopo_griddedInterpolant, ...
        longwaveTopo_griddedInterpolant,Topo_PARA.Density,Coastline)
else
    disp('2/4 ..........................computeTerrainEffect is running')
    [fullTopoCorrectedGravityPoint, longwaveTopo_griddedInterpolant, fullTopo_griddedInterpolant, fullTopoCorrectedGravityGradient] = ...
        computeFullTerrainEffects(GRID_PARA, Topo_PARA, Gravo, gravGradFiltered, GGM_Gravity_griddedInterpolant, DEM_data, ZDEM_griddedInterpolant, ...
        LongDEM, LatDEM, Coastline, OUTPUT_PARA.plotsFolder);

    save([OUTPUT_PARA.Grids_name, 'terrainEffects','.mat'], 'fullTopoCorrectedGravityPoint', 'longwaveTopo_griddedInterpolant', 'fullTopo_griddedInterpolant', 'fullTopoCorrectedGravityGradient');

    disp('3/4 ..........................computeGravimetryGradiometryLSC is running')
    computeGravimetryGradiometryLSC(GRID_PARA,COV_PARA,DEM_PARA,GRAV_PARA,GRAV_GRAD_PARA,OUTPUT_PARA,GRID_REF,fullTopoCorrectedGravityPoint,fullTopoCorrectedGravityGradient, ...
        GGM_Gravity_griddedInterpolant,ZDEM_griddedInterpolant,fullTopo_griddedInterpolant, ...
        longwaveTopo_griddedInterpolant,Topo_PARA.Density,Coastline)
end

disp('4/4 ..........................mosaicTiles is running')
geomGravGeoidDiff = mosaicTiles(GRID_PARA,DEM_PARA,OUTPUT_PARA,Lev,LongDEM,LatDEM, ...
    REFERENCE_Zeta_griddedInterpolant,GGM_Gravity_griddedInterpolant,GGM_Zeta_griddedInterpolant,Coastline);

diary off
