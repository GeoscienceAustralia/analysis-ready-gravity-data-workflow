% Gravimetric quasigeoid estimation
%
% RunDemo is the main m-file preforming the gravimetric quasiGeoid estimation.
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
%             - computeGravimetryLSC: runs the LSC (Least Squares Collocation) in blocks.
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
GRID_PARA.buffer=1;% degs. The x/y extent to extract data around the tile centre. .75
GRID_PARA.buffer2=0.5;% degs. The x/y tile extents that are kept - where the good data are.
GRID_PARA.STEP=0.5;% The step size. This must be less than buffer2 to avoid gaps in the final grid.
GRID_PARA.filterSize=15;% filter size for spatial grid weight, this value is from experiment for tiles of one degree
GRID_PARA.filterRadius=10; % filter radius for spatial grid weight, this value is from experiment for tiles of one degree
% Grid extents - ensure these values are in GRID_PARA.STEP degree value increments.
% Boundary for computation
GRID_PARA.MINLONG=141.5;
GRID_PARA.MAXLONG=142;
GRID_PARA.MINLAT=-32.5;
GRID_PARA.MAXLAT=-32;
%% DEM data - N.B. the dem is used to specify the grid nodes.
DEM_PARA.filename='Data\DEM\AUSDEM1min.xyz';
DEM_PARA.num_cols=4861; % original size, to be changed, see below
DEM_PARA.num_rows=3181; % original size, to be changed, see below
%% Gravity data
% this collates all of the gravity and position data into one matlab array.
GRAV_PARA.filename='Data\processedData\Terrestrial_Gravity.mat';
GRAV_PARA.filename1=[];%'Data\GRAVITY\Xcalibur_Gravity.mat';
GRAV_PARA.TypeB=1;% This is a Type B uncertainty value (in mGal) which is added to the uncertainty values.
GRAV_PARA.Grav_Faye_TypeB=3;
%% Gravity Gradiometry data
GRAV_GRAD_PARA.filename='Data\GRAVITY_GRAD\Xcalibur_FVD_GDD.mat';
GRAV_GRAD_PARA.TypeB=10^(-5);% This is a Type B uncertainty value (in mGal/m) which is added to the uncertainty values.
GRAV_GRAD_PARA.avail=false;
%% Covariance function parameters
COV_PARA.Compute_Empircal_COV_Dec=3; % Decimation factor for empirical covariance estimates. e.g. 1 is no decimation, 2 drops 50% of the data etc. see sph_empcov for logic.
COV_PARA.Fit_Empircal_COV='auto';%'auto';% process to fit covariance N & M function values 'man' for manual to fit them on the cmd line or 'auto' , '' to just use what you supply here.
COV_PARA.FitEmpiricalCOVNSearch=[21600,1,21600];% Start, step, stop parameter sweep values for N parameter - if auto
COV_PARA.FitEmpiricalCOVMSearch=[200,20,300];% Start, step, stop parameter sweep values for M parameter - if auto
COV_PARA.N=10800;% max Legender polynonial of cov func. 
COV_PARA.M=200;% min Legender polynonial of cov func. 
COV_PARA.width=3;% Size of precomputed cov function in degrees - must be larger the the distance between any two points on a tile. 
COV_PARA.res=30/3600; % Resolution of the covariance function
COV_PARA.COV_COMPUTED_Tilewise=true;% This recomputes the covariance function for each tile.
COV_PARA.Airbornedataonly=false;% Only use airborne data in establishing Covariance parameters - good to use if we are using EGM2008 as the references as terrestrial data are not independent.
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
GGM_PARA.filename='Data\processedData\GOCE_N300_For_Gridded_Int24_2height.mat';
%% Coastline data
COAST_PARA.filename='Data\COASTLINE\CoastAus.mat';
%% Levelling data comparisons
LEVELLING_PARA.Lev_eval=true;% If true, the levelling data are compared to the geoid as its computed.
LEVELLING_PARA.filename='Data\GPS_LEVELLING\Lev_NSW_NG.mat';%'Data/GPS_LEVELLING/Lev_CARS.mat';% The format of these data needs to be an array with rows [Long,Lat,h-H].
LEVELLING_PARA.Plot_Stats=true;% If true, the levelling data are compared to the geoid as its computed.
LEVELLING_PARA.Compare_To_Existing_Model=true;% If true, the levelling data are also compared to another existing geoid as its computed.
LEVELLING_PARA.Existing_Model='Data\EXISTING_GEOID_MODELS\AGQG20221120.mat';% File location of the existing model.
LEVELLING_PARA.max_diff=0.15;% Threshold for an outlier with the GNSS-levelling
%% Output
OUTPUT_PARA.Grids_name = 'outputsDemo/Grids/';
OUTPUT_PARA.Tiles_dir_name = 'outputsDemo/ResidualTiles/';
OUTPUT_PARA.PLOT_GRIDS = true;% A gridded solution is plotted and output as well as the tiles.
OUTPUT_PARA.plotsFolder = 'outputsDemo/plots/';
% If there is a region of interest, for plotting purposes
OUTPUT_PARA.polygonLon = [];
OUTPUT_PARA.polygonLat = [];

% Keep the computer awake
keepawake=true;% Setting this to true wiggles the mouse every so often so the compute doesnt go to sleep.
%% Run the LSC code
% this is how demo data is created
% disp('1/4 ..........................importAndFormatData is running ')
% [Gravity6D,gravGradFiltered,DEM3D,ZDEM_griddedInterpolant,LongDEMmatrix,LatDEMmatrix,...
%  GravityGGM_griddedInterpolant,ZetaGGM_griddedInterpolant,LevellingData3D,...
%  ZetaRef_griddedInterpolant,GridRef3D,Coastline,DEM_PARA]=importAndFormatData...
%  (GRID_PARA,DEM_PARA,GRAV_PARA,Topo_PARA,COAST_PARA,LEVELLING_PARA,GGM_PARA,GRAV_GRAD_PARA);
% 
% save([OUTPUT_PARA.Grids_name,'oneTileData141-32','.mat'],...
%  'Gravity6D','DEM3D','LongDEMmatrix','LatDEMmatrix','LevellingData3D','GridRef3D','Coastline')
% 
% save([OUTPUT_PARA.Grids_name,'oneTileGriddedInterpolant141-32','.mat'],...
% 'ZDEM_griddedInterpolant','GravityGGM_griddedInterpolant','ZetaGGM_griddedInterpolant','ZetaRef_griddedInterpolant')

% Check if the directory exists, if not, create it
if ~exist(OUTPUT_PARA.Grids_name, 'dir')
    mkdir(OUTPUT_PARA.Grids_name);
end

if ~exist(OUTPUT_PARA.Tiles_dir_name, 'dir')
    mkdir(OUTPUT_PARA.Tiles_dir_name);
end

if ~exist(OUTPUT_PARA.plotsFolder, 'dir')
    mkdir(OUTPUT_PARA.plotsFolder);
end

disp('1/4 ..........................import demo data')

dataDemo=importdata('Data\oneTileData141-32.mat');
interpolantDemo=importdata('Data\oneTileGriddedInterpolant141-32.mat');

DEM_PARA.num_rows=size(dataDemo.LatDEMmatrix,1);
DEM_PARA.num_cols=size(dataDemo.LatDEMmatrix,2);

% Plot input data: 

gravGradFiltered=[];

if OUTPUT_PARA.PLOT_GRIDS
     plotInputData(dataDemo.Gravity6D,gravGradFiltered,dataDemo.Coastline,GRID_PARA,OUTPUT_PARA)
end 


if exist([OUTPUT_PARA.Grids_name,'terrainEffects.mat'], 'file')
    load([OUTPUT_PARA.Grids_name,'terrainEffects.mat']);
    disp('3/4 ..........................computeGravimetryLSC is running')
    computeGravimetryLSC(GRID_PARA,COV_PARA,DEM_PARA,GRAV_PARA,OUTPUT_PARA,dataDemo.GridRef3D,fullTopoCorrectedGravityPoint, ...
    interpolantDemo.GravityGGM_griddedInterpolant,interpolantDemo.ZDEM_griddedInterpolant,fullTopo_griddedInterpolant, ...
    longwaveTopo_griddedInterpolant,Topo_PARA.Density)

else
    disp('2/4 ..........................computeTerrainEffect is running')
    [fullTopoCorrectedGravityPoint,longwaveTopo_griddedInterpolant,fullTopo_griddedInterpolant]=computeTerrainEffect(GRID_PARA, ...
    Topo_PARA,dataDemo.Gravity6D,interpolantDemo.GravityGGM_griddedInterpolant,dataDemo.DEM3D,interpolantDemo.ZDEM_griddedInterpolant, ...
    dataDemo.LongDEMmatrix,dataDemo.LatDEMmatrix,dataDemo.Coastline,OUTPUT_PARA.plotsFolder);

    save([OUTPUT_PARA.Grids_name,'TerrainEffects','.mat'],'fullTopoCorrectedGravityPoint','longwaveTopo_griddedInterpolant','fullTopo_griddedInterpolant')

    disp('3/4 ..........................computeGravimetryLSC is running')
    computeGravimetryLSC(GRID_PARA,COV_PARA,DEM_PARA,GRAV_PARA,OUTPUT_PARA,dataDemo.GridRef3D,fullTopoCorrectedGravityPoint, ...
    interpolantDemo.GravityGGM_griddedInterpolant,interpolantDemo.ZDEM_griddedInterpolant,fullTopo_griddedInterpolant, ...
    longwaveTopo_griddedInterpolant,Topo_PARA.Density)
end

disp('4/4 ..........................mosaicTiles is running')
geomGravGeoidDiff = mosaicTiles(GRID_PARA,DEM_PARA,OUTPUT_PARA,dataDemo.LevellingData3D,dataDemo.LongDEMmatrix,dataDemo.LatDEMmatrix, ...
    interpolantDemo.ZetaRef_griddedInterpolant,interpolantDemo.GravityGGM_griddedInterpolant,interpolantDemo.ZetaGGM_griddedInterpolant,dataDemo.Coastline);

% % test Parallel pull for tilewise LSC 
% tic
% 
% disp('3/4 ..........................computeParallelGravimetryLSC')
% computeParallelGravimetryLSC(GRID_PARA,COV_PARA,DEM_PARA,GRAV_PARA,OUTPUT_PARA,dataDemo.GridRef3D,fullTopoCorrectedGravityPoint, ...
%     interpolantDemo.GravityGGM_griddedInterpolant,interpolantDemo.ZDEM_griddedInterpolant,fullTopo_griddedInterpolant, ...
%     longwaveTopo_griddedInterpolant,Topo_PARA.Density)
% toc

