function RunParallelMainScript(varargin)
%RunMainScript computes regional gravimetric geoids using gravity observations from gravity anomalies.
%              The process involves sequence of "remove-predict-restore" operations, where the Global
%              Gravity Model (GGM) and topographic effects are removed, a geoid is predicted (here with LSC),
%              and then the effects are restored to obtain the final geoid model. The functions folder
%              provides all the MATLAB functions to perform these three steps for geoid calculations.
%              The primary goal is to create a platform for analysis-ready gravity data,
%              featuring a tile-wise least-squares collocation (LSC) method based on gravity anomaly
%              observations.
%
% Usage: RunMainScript('flag',value)
%    or: RunMainScript flag value
%
%e.g.
%Within Matlab
%             RunMainScript('--grid-para-buffer',1,'--dem-para-filename','/g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz,'--grav-grad-para-avail',true);
%
%Compiled
%             RunMainScript --grid-para-buffer 1 --dem-para-filename /g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz --grav-grad-para-avail true
%
%Available options:
%--grid-para-buffer <value>                             e.g. --grid-para-buffer 1
%--grid-para-buffer2 <value>                            e.g. --grid-para-buffer2 0.5
%--grid_para-step <value>                               e.g. --grid_para-step 0.5
%--grid-para-filtersize <value>                         e.g. --grid-para-filtersize 15
%--grid-para-filterradius <value>                       e.g. --grid-para-filterradius 10
%--grid-para-minlong <value>                            e.g. --grid-para-minlong 153
%--grid-para-maxlong <value>                            e.g. --grid-para-maxlong 154
%--grid-para-minlat <value>                             e.g. --grid-para-minlat -29
%--grid-para-maxlat <value>                             e.g. --grid-para-maxlat -28
%--dem-para-filename <path_to_file>                     e.g. --dem-para-filename /g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz
%--dem-para-num-cols <value>                            e.g. --dem-para-num-cols 4861
%--dem-para-num-rows <value>                            e.g. --dem-para-num-rows 3181
%--grav-para-filename <path_to_file>                    e.g. --grav-para-filename /g/data/dg9/nd2979/Data/processedData/GravityAllVicNSW.mat
%--grav-para-filename1 <path_to_file>                   e.g. --grav-para-filename1 /g/data/dg9/nd2979/Data/processedData/GravityAllVicNSW_1.mat
%--grav-para-typeb <value>                              e.g. --grav-para-typeb 1
%--grav-para-grav-faye-typeb <value>                    e.g. --grav-para-grav-faye -typeb 3
%--grav-para-altimetry-weighting <value>                e.g. --grav-para-altimetry-weighting 1
%--grav-grad-para-filename <path_to_file>               e.g. --grav-grad-para-filename /g/data/dg9/nd2979/Data/GRAVITY_GRAD/Xcalibur_FVD_GDD.mat
%--grav-grad-para-typeb <value>                         e.g. --grav-grad-para-typeb 0.00001
%--grav-grad-para-avail <logical>                       e.g. --grav-grad-para-avail true
%--cov-para-compute-empircal-cov-dec <value>            e.g. --cov-para-compute-empircal-cov-dec 3
%--cov-para-fit-empircal-cov <type>                     e.g. --cov-para-fit-empircal-cov auto
%--cov-para-fitempiricalcovnsearch <values>             e.g. --cov-para-fitempiricalcovnsearch 21600,1,21600
%--cov-para-fitempiricalcovmsearch <values>             e.g. --cov-para-fitempiricalcovmsearch 200,20,300
%--cov-para-n <value>                                   e.g. --cov-para-n 10800
%--cov-para-m <value>                                   e.g. --cov-para-m 200
%--cov-para-width <value>                               e.g. --cov-para-width 3
%--cov-para-res <value>                                 e.g. --cov-para-res 0.00833333333
%--cov-para-cov-computed_tilewise <logical>             e.g. --cov-para-cov-computed_tilewise true
%--cov-para-airbornedataonly <logical>                  e.g. --cov-para-airbornedataonly false
%--cov-para-covplot <logical>                           e.g. --cov-para-covplot false
%--topo-para-corr <logical>                             e.g. --topo-para-corr true
%--topo-para-topoplot <logical>                         e.g. --topo-para-topoplot false
%--topo-para-density <value>                            e.g. --topo-para-density 2.67
%--topo-para-depth <value>                              e.g. --topo-para-depth 0
%--topo-para-rad <value>                                e.g. --topo-para-rad 1
%--topo-para-rtm <values>                               e.g. --topo-para-rtm 50,10,300
%--ggm-para-filename <path_to_file>                     e.g. --ggm-para-filename /g/data/dg9/nd2979/Data/GGM/GOCE_For_Gridded_Int.mat
%--coast-para-filename <path_to_file>                   e.g. --coast-para-filename /g/data/dg9/nd2979/Data/COASTLINE/CoastAus.mat
%--levelling-para-lev-eval <logical>                    e.g. --levelling-para-lev-eval true
%--levelling-para-filename <path_to_file>               e.g. --levelling-para-filename /g/data/dg9/nd2979/Data/GPS_LEVELLING/Lev_NSW_NG.mat
%--levelling-para-plot-stats <logical>                  e.g. --levelling-para-plot-stats false
%--levelling-para-compare-to-existing-model <logical>   e.g. --levelling-para-compare-to-existing-model true
%--levelling-para-existing-model <path_to_file>         e.g. --levelling-para-existing-model /g/data/dg9/nd2979/Data/EXISTING_GEOID_MODELS/AGQG20221120.mat
%--levelling-para-max-diff <value>                      e.g. --levelling-para-max-diff 0.15
%--output-para-grids-name <path_to_folder>              e.g. --output-para-grids-name /g/data/dg9/nd2979/outputs/GridsNENSWgg2degTile/
%--output-para-tiles-dir-name <path_to_folder>          e.g. --output-para-tiles-dir-name /g/data/dh8/outputs/ResidualTilesNENSWgg2degTile/
%--output-para-plot-grids <logical>                     e.g. --output-para-plot-grids false
%--output-para-plotsfolder <path_to_folder>             e.g. --output-para-plotsfolder /g/data/dh8/outputs/plots/22-Nov-2024NENSWgg2degTile
%--keepawake <logical>                                  e.g. --keepawake true
%
%Geoscience Australia. Neda Darbeheshti on 24/11/2024
%

% Grid/Tiling Parameters
% Tiling Parameters - fixed for each computation
GRID_PARA.buffer=1; % degs. The x/y extent to extract data around the tile centre. .75
GRID_PARA.buffer2=0.5; % degs. The x/y tile extents that are kept - where the good data are.
GRID_PARA.STEP=0.5; % The step size. This must be less than buffer2 to avoid gaps in the final grid.
GRID_PARA.filterSize=15; % filter size for spatial grid weight, this value is from experiment for tiles of one degree
GRID_PARA.filterRadius=10; % filter radius for spatial grid weight, this value is from experiment for tiles of one degree
% Grid extents - ensure these values are in GRID_PARA.STEP degree value increments.
% Boundary for computation
GRID_PARA.MINLONG=153;
GRID_PARA.MAXLONG=154;
GRID_PARA.MINLAT=-29;
GRID_PARA.MAXLAT=-28;

% DEM data - N.B. the dem is used to specify the grid nodes.
DEM_PARA.filename='/g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz';
DEM_PARA.num_cols=4861; % original size, to be changed, see below
DEM_PARA.num_rows=3181; % original size, to be changed, see below

% Gravity data
GRAV_PARA.filename='/g/data/dg9/nd2979/Data/processedData/GravityAllVicNSW.mat';
GRAV_PARA.filename1=[]; %'Data\GRAVITY\Xcalibur_Gravity.mat';
GRAV_PARA.TypeB=1; % This is a Type B uncertainty value (in mGal) which is added to the uncertainty values.
GRAV_PARA.Grav_Faye_TypeB=3;
GRAV_PARA.inputGravity_weighting=0;

% Gravity Gradiometry data
GRAV_GRAD_PARA.filename='/g/data/dg9/nd2979/Data/GRAVITY_GRAD/Xcalibur_FVD_GDD.mat';
GRAV_GRAD_PARA.TypeB=10^(-5); % This is a Type B uncertainty value (in mGal/m) which is added to the uncertainty values.
GRAV_GRAD_PARA.avail=true;

% Covariance function parameters
COV_PARA.Compute_Empircal_COV_Dec=3; % Decimation factor for empirical covariance estimates. e.g. 1 is no decimation, 2 drops 50% of the data etc. see sph_empcov for logic.
COV_PARA.Fit_Empircal_COV='auto'; % process to fit covariance N & M function values 'man' for manual to fit them on the cmd line or 'auto' , '' to just use what you supply here.
COV_PARA.FitEmpiricalCOVNSearch=[21600,1,21600]; % Start, step, stop parameter sweep values for N parameter - if auto
COV_PARA.FitEmpiricalCOVMSearch=[200,20,300]; % Start, step, stop parameter sweep values for M parameter - if auto
COV_PARA.N=10800; % max Legender polynonial of cov func. 
COV_PARA.M=200; % min Legender polynonial of cov func. 
COV_PARA.width=3; % Size of precomputed cov function in degrees - must be larger the the distance between any two points on a tile. 
COV_PARA.res=30/3600; % Resolution of the covariance function
COV_PARA.COV_COMPUTED_Tilewise=true; % This recomputes the covariance function for each tile.
COV_PARA.Airbornedataonly=false; % Only use airborne data in establishing Covariance parameters - good to use if we are using EGM2008 as the references as terrestrial data are not independent.
COV_PARA.COVPlot=false; % true plots progress, false turns this off.

% Topo condensation parameters
Topo_PARA.Corr=true; % MAKE SURE YOU TURN THIS ON!!!
Topo_PARA.TopoPlot=false; % true plots progress, false turns this off.
Topo_PARA.Density=2.67; % Assumed density in g/cm^3.
Topo_PARA.Depth=0; % Condensation layer depth. 0 is on the geoid
Topo_PARA.Rad=1; % Radius out to which to compute the effects in degress.
Topo_PARA.RTM=[50,10,300]; %[1000,10,2160] for egm%[0,10,300] % Range of SHM degree filter parameters (min, step, max) explored when running RTM calculations.

% GGM reference signal
GGM_PARA.filename='/g/data/dg9/nd2979/Data/GGM/GOCE_For_Gridded_Int.mat';

% Coastline data
COAST_PARA.filename='/g/data/dg9/nd2979/Data/COASTLINE/CoastAus.mat';

% Levelling data comparisons
LEVELLING_PARA.Lev_eval=true; % If true, the levelling data are compared to the geoid as its computed.
LEVELLING_PARA.filename='/g/data/dg9/nd2979/Data/GPS_LEVELLING/Lev_NSW_NG.mat'; %'/g/data/dg9/nd2979/Data/GPS_LEVELLING/Lev_CARS.mat';% The format of these data needs to be an array with rows [Long,Lat,h-H].
LEVELLING_PARA.Plot_Stats=false; % If true, the levelling data are compared to the geoid as its computed.
LEVELLING_PARA.Compare_To_Existing_Model=true; % If true, the levelling data are also compared to another existing geoid as its computed.
LEVELLING_PARA.Existing_Model='/g/data/dg9/nd2979/Data/EXISTING_GEOID_MODELS/AGQG20221120.mat'; % File location of the existing model.
LEVELLING_PARA.max_diff=0.15;% Threshold for an outlier with the GNSS-levelling

% Output
OUTPUT_PARA.Grids_name='/g/data/dg9/nd2979/outputs/GridsNENSWgg2degTile/';
OUTPUT_PARA.Tiles_dir_name='/g/data/dg9/nd2979/outputs/ResidualTilesNENSWgg2degTile/';
OUTPUT_PARA.PLOT_GRIDS=false;% A gridded solution is plotted and output as well as the tiles.
OUTPUT_PARA.plotsFolder='/g/data/dg9/nd2979/outputs/plots/22-Nov-2024NENSWgg2degTile';
% If there is a region of interest, for plotting purposes
OUTPUT_PARA.polygonLon = [115.4333, 116.0500, 116.2500, 116.2500, 115.6167, 115.6167, 115.4333 ];
OUTPUT_PARA.polygonLat = [-31.4500, -31.4500, -32.0000, -32.5833, -32.5833, -32.0000,-31.4500];

% Keep the computer awake
keepawake=true; % Setting this to true wiggles the mouse every so often so the compute doesnt go to sleep.

%check for input argument
for i=1:nargin
    if ischar(varargin{i})
        if strncmp(varargin{i},'--grid-para-buffer',18)
            GRID_PARA.buffer = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grid-para-buffer2',19)
            GRID_PARA.buffer2 = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grid_para-step',16)
            GRID_PARA.STEP = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grid-para-filtersize',22)
            GRID_PARA.filterSize = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grid-para-filterradius',24)
            GRID_PARA.filterRadius = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grid-para-minlong',19)
            GRID_PARA.MINLONG = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grid-para-maxlong',19)
            GRID_PARA.MAXLONG = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grid-para-minlat',18)
            GRID_PARA.MINLAT = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grid-para-maxlat',18)
            GRID_PARA.MAXLAT = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--dem-para-filename',19)
            DEM_PARA.filename = varargin{i+1};
        elseif strncmp(varargin{i},'--dem-para-num-cols',19)
            DEM_PARA.num_cols = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--dem-para-num-rows',19)
            DEM_PARA.num_rows = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grav-para-filename',20)
            GRAV_PARA.filename = varargin{i+1};
        elseif strncmp(varargin{i},'--grav-para-filename1',21)
            GRAV_PARA.filename1 = varargin{i+1};
        elseif strncmp(varargin{i},'--grav-para-typeb',17)
            GRAV_PARA.TypeB = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grav-para-grav-faye-typeb',27)
            GRAV_PARA.Grav_Faye_TypeB = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grav-para-altimetry-weighting',21)
            GRAV_PARA.altimetry_weighting = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grav-grad-para-filename',25)
            GRAV_GRAD_PARA.filename = varargin{i+1};
        elseif strncmp(varargin{i},'--grav-grad-para-typeb',22)
            GRAV_GRAD_PARA.TypeB = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--grav-grad-para-avail',22)
           GRAV_GRAD_PARA.avail = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--cov-para-compute-empircal-cov-dec',35)
           COV_PARA.Compute_Empircal_COV_Dec = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--cov-para-fit-empircal-cov',27)
           COV_PARA.Fit_Empircal_COV = varargin{i+1};
        elseif strncmp(varargin{i},'--cov-para-fitempiricalcovnsearch',33)
           COV_PARA.FitEmpiricalCOVNSearch = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--cov-para-fitempiricalcovmsearch',33)
           COV_PARA.FitEmpiricalCOVMSearch = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--cov-para-n',12)
           COV_PARA.N = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--cov-para-m',12)
           COV_PARA.M = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--cov-para-width',16)
          COV_PARA.width = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--cov-para-res',14)
          COV_PARA.res = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--cov-para-cov-computed_tilewise',32)
          COV_PARA.COV_COMPUTED_Tilewise = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--cov-para-airbornedataonly',27)
          COV_PARA.Airbornedataonly = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--cov-para-covplot',18)
          COV_PARA.COVPlot = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--topo-para-corr',16)
          Topo_PARA.Corr = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--topo-para-topoplot',20)
          Topo_PARA.TopoPlot = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--topo-para-density',19)
          Topo_PARA.Density = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--topo-para-depth',17)
          Topo_PARA.Depth = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--topo-para-rad',15)
          Topo_PARA.Rad = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--topo-para-rtm',15)
          Topo_PARA.RTM = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--ggm-para-filename',19)
          GGM_PARA.filename = varargin{i+1};
        elseif strncmp(varargin{i},'--coast-para-filename',21)
          COAST_PARA.filename = varargin{i+1};
        elseif strncmp(varargin{i},'--levelling-para-lev-eval',25)
          LEVELLING_PARA.Lev_eval = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--levelling-para-filename',25)
          LEVELLING_PARA.filename = varargin{i+1};
        elseif strncmp(varargin{i},'--levelling-para-plot-stats',27)
          LEVELLING_PARA.Plot_Stats = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--levelling-para-compare-to-existing-model',42)
          LEVELLING_PARA.Compare_To_Existing_Model = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--levelling-para-existing-model',31)
          LEVELLING_PARA.Existing_Model = varargin{i+1};
        elseif strncmp(varargin{i},'--levelling-para-max-diff',25)
          LEVELLING_PARA.max_diff = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--output-para-grids-name',24)
          OUTPUT_PARA.Grids_name = varargin{i+1};
        elseif strncmp(varargin{i},'--output-para-tiles-dir-name',28)
          OUTPUT_PARA.Tiles_dir_name = varargin{i+1};
        elseif strncmp(varargin{i},'--output-para-plot-grids',24)
          OUTPUT_PARA.PLOT_GRIDS = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--output-para-plotsfolder',25)
          OUTPUT_PARA.plotsFolder = varargin{i+1};
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strncmp(varargin{i},'--output-para-polygonLon',24)
          OUTPUT_PARA.polygonLon = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--output-para-polygonLat',24)
          OUTPUT_PARA.polygonLat = str2num(varargin{i+1});
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strncmp(varargin{i},'--keepawake',11)
          keepawake = str2num(varargin{i+1});
        elseif strncmp(varargin{i},'--help',6)
          helptext
          return
        end
    end
end

%check for the presence of input file
if ~isfile(DEM_PARA.filename)
    disp(['File ' DEM_PARA.filename ' not found'])
    return
end
if ~isfile(GRAV_PARA.filename)
    disp(['File ' GRAV_PARA.filename ' not found'])
    return
end
if ~isfile(GRAV_GRAD_PARA.filename)
    disp(['File ' GRAV_GRAD_PARA.filename ' not found'])
    return
end
if ~isfile(GGM_PARA.filename)
    disp(['File ' GGM_PARA.filename ' not found'])
    return
end
if ~isfile(COAST_PARA.filename)
    disp(['File ' COAST_PARA.filename ' not found'])
    return
end
if ~isfile(LEVELLING_PARA.filename)
    disp(['File ' LEVELLING_PARA.filename ' not found'])
    return
end
if ~isfile(LEVELLING_PARA.Existing_Model)
    disp(['File ' LEVELLING_PARA.Existing_Model ' not found'])
    return
end

%check for output folders, create as necessary
if ~isfolder(OUTPUT_PARA.Grids_name)
    mkdir(OUTPUT_PARA.Grids_name)
end
if ~isfolder(OUTPUT_PARA.Tiles_dir_name)
    mkdir(OUTPUT_PARA.Tiles_dir_name)
end
if ~isfolder(OUTPUT_PARA.plotsFolder)
    mkdir(OUTPUT_PARA.plotsFolder)
end

%temporary breakpoint to check
%disp(GRID_PARA)
%disp(GRAV_GRAD_PARA)
%disp(OUTPUT_PARA)
%return

%computation
disp('1/4 ..........................importAndFormatData is running ')
[Gravo,gravGradFiltered,DEM_data,ZDEM_griddedInterpolant,LongDEM,LatDEM,...
 GGM_Gravity_griddedInterpolant,GGM_Zeta_griddedInterpolant,Lev,...
 REFERENCE_Zeta_griddedInterpolant,GRID_REF,Coastline,DEM_PARA]=importAndFormatData...
 (GRID_PARA,DEM_PARA,GRAV_PARA,Topo_PARA,COAST_PARA,LEVELLING_PARA,GGM_PARA,GRAV_GRAD_PARA);

if OUTPUT_PARA.PLOT_GRIDS
     plotInputData(Gravo,gravGradFiltered,Coastline,GRID_PARA,OUTPUT_PARA)
end 

if GRAV_PARA.inputGravity_weighting 
     Gravo = weightInputGravity(Gravo,Coastline,GRID_PARA,OUTPUT_PARA);
end

if exist([OUTPUT_PARA.Grids_name,'terrainEffects.mat'], 'file')
    load([OUTPUT_PARA.Grids_name,'terrainEffects.mat']);
    disp('3/4 ..........................computeGravimetryGradiometryLSC is running')
    computeParallelGravimetryGradiometryLSC(GRID_PARA,COV_PARA,DEM_PARA,GRAV_PARA,GRAV_GRAD_PARA,OUTPUT_PARA,GRID_REF,fullTopoCorrectedGravityPoint,fullTopoCorrectedGravityGradient, ...
        GGM_Gravity_griddedInterpolant,ZDEM_griddedInterpolant,fullTopo_griddedInterpolant, ...
        longwaveTopo_griddedInterpolant,Topo_PARA.Density,Coastline)

else
    disp('2/4 ..........................computeTerrainEffect is running')
    [fullTopoCorrectedGravityPoint, longwaveTopo_griddedInterpolant, fullTopo_griddedInterpolant, fullTopoCorrectedGravityGradient] = ...
        computeFullTerrainEffects(GRID_PARA, Topo_PARA, Gravo, gravGradFiltered, GGM_Gravity_griddedInterpolant, DEM_data, ZDEM_griddedInterpolant, ...
        LongDEM, LatDEM, Coastline, OUTPUT_PARA.plotsFolder);

    save([OUTPUT_PARA.Grids_name, 'terrainEffects','.mat'], 'fullTopoCorrectedGravityPoint', 'longwaveTopo_griddedInterpolant', 'fullTopo_griddedInterpolant', 'fullTopoCorrectedGravityGradient');

    disp('3/4 ..........................computeGravimetryGradiometryLSC is running')
    computeParallelGravimetryGradiometryLSC(GRID_PARA,COV_PARA,DEM_PARA,GRAV_PARA,GRAV_GRAD_PARA,OUTPUT_PARA,GRID_REF,fullTopoCorrectedGravityPoint,fullTopoCorrectedGravityGradient, ...
        GGM_Gravity_griddedInterpolant,ZDEM_griddedInterpolant,fullTopo_griddedInterpolant, ...
        longwaveTopo_griddedInterpolant,Topo_PARA.Density,Coastline)
end

% disp('4/4 ..........................mosaicTiles is running')
% geomGravGeoidDiff = mosaicTiles(GRID_PARA,DEM_PARA,OUTPUT_PARA,Lev,LongDEM,LatDEM, ...
%     REFERENCE_Zeta_griddedInterpolant,GGM_Gravity_griddedInterpolant,GGM_Zeta_griddedInterpolant,Coastline);

function helptext
str={'RunMainScript computes regional gravimetric geoids using gravity observations from gravity anomalies.'
'              The process involves sequence of "remove-predict-restore" operations, where the Global '
'              Gravity Model (GGM) and topographic effects are removed, a geoid is predicted (here with LSC), '
'              and then the effects are restored to obtain the final geoid model. The functions folder '
'              provides all the MATLAB functions to perform these three steps for geoid calculations. '
'              The primary goal is to create a platform for analysis-ready gravity data, '
'              featuring a tile-wise least-squares collocation (LSC) method based on gravity anomaly '
'              observations.'
' '
' Usage: RunMainScript(''flag'',value)'
'    or: RunMainScript flag value'
' '
'e.g.'
'Within Matlab'
'             RunMainScript(''--grid-para-buffer'',1,''--dem-para-filename'',''/g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz,''--grav-grad-para-avail'',true);'
' '
'Compiled'
'             RunMainScript --grid-para-buffer 1 --dem-para-filename /g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz --grav-grad-para-avail true'
' '
'Available options:'
'--grid-para-buffer <value>                             e.g. --grid-para-buffer 1'
'--grid-para-buffer2 <value>                            e.g. --grid-para-buffer2 0.5'
'--grid_para-step <value>                               e.g. --grid_para-step 0.5'
'--grid-para-filtersize <value>                         e.g. --grid-para-filtersize 15'
'--grid-para-filterradius <value>                       e.g. --grid-para-filterradius 10'
'--grid-para-minlong <value>                            e.g. --grid-para-minlong 153'
'--grid-para-maxlong <value>                            e.g. --grid-para-maxlong 154'
'--grid-para-minlat <value>                             e.g. --grid-para-minlat -29'
'--grid-para-maxlat <value>                             e.g. --grid-para-maxlat -28'
'--dem-para-filename <path_to_file>                     e.g. --dem-para-filename /g/data/dg9/nd2979/Data/DEM/AUSDEM1min.xyz'
'--dem-para-num-cols <value>                            e.g. --dem-para-num-cols 4861'
'--dem-para-num-rows <value>                            e.g. --dem-para-num-rows 3181'
'--grav-para-filename <path_to_file>                    e.g. --grav-para-filename /g/data/dg9/nd2979/Data/processedData/GravityAllVicNSW.mat'
'--grav-para-filename1 <path_to_file>                   e.g. --grav-para-filename1 /g/data/dg9/nd2979/Data/processedData/GravityAllVicNSW_1.mat'
'--grav-para-typeb <value>                              e.g. --grav-para-typeb 1'
'--grav-para-grav-faye-typeb <value>                    e.g. --grav-para-grav-faye -typeb 3'
'--grav-para-altimetry-weighting <value>                e.g. --grav-para-altimetry-weighting 1'
'--grav-grad-para-filename <path_to_file>               e.g. --grav-grad-para-filename /g/data/dg9/nd2979/Data/GRAVITY_GRAD/Xcalibur_FVD_GDD.mat'
'--grav-grad-para-typeb <value>                         e.g. --grav-grad-para-typeb 0.00001'
'--grav-grad-para-avail <logical>                       e.g. --grav-grad-para-avail true'
'--cov-para-compute-empircal-cov-dec <value>            e.g. --cov-para-compute-empircal-cov-dec 3'
'--cov-para-fit-empircal-cov <type>                     e.g. --cov-para-fit-empircal-cov auto'
'--cov-para-fitempiricalcovnsearch <values>             e.g. --cov-para-fitempiricalcovnsearch 21600,1,21600'
'--cov-para-fitempiricalcovmsearch <values>             e.g. --cov-para-fitempiricalcovmsearch 200,20,300'
'--cov-para-n <value>                                   e.g. --cov-para-n 10800'
'--cov-para-m <value>                                   e.g. --cov-para-m 200'
'--cov-para-width <value>                               e.g. --cov-para-width 3'
'--cov-para-res <value>                                 e.g. --cov-para-res 0.00833333333'
'--cov-para-cov-computed_tilewise <logical>             e.g. --cov-para-cov-computed_tilewise true'
'--cov-para-airbornedataonly <logical>                  e.g. --cov-para-airbornedataonly false'
'--cov-para-covplot <logical>                           e.g. --cov-para-covplot false'
'--topo-para-corr <logical>                             e.g. --topo-para-corr true'
'--topo-para-topoplot <logical>                         e.g. --topo-para-topoplot false'
'--topo-para-density <value>                            e.g. --topo-para-density 2.67'
'--topo-para-depth <value>                              e.g. --topo-para-depth 0'
'--topo-para-rad <value>                                e.g. --topo-para-rad 1'
'--topo-para-rtm <values>                               e.g. --topo-para-rtm 50,10,300'
'--ggm-para-filename <path_to_file>                     e.g. --ggm-para-filename /g/data/dg9/nd2979/Data/GGM/GOCE_For_Gridded_Int.mat'
'--coast-para-filename <path_to_file>                   e.g. --coast-para-filename /g/data/dg9/nd2979/Data/COASTLINE/CoastAus.mat'
'--levelling-para-lev-eval <logical>                    e.g. --levelling-para-lev-eval true'
'--levelling-para-filename <path_to_file>               e.g. --levelling-para-filename /g/data/dg9/nd2979/Data/GPS_LEVELLING/Lev_NSW_NG.mat'
'--levelling-para-plot-stats <logical>                  e.g. --levelling-para-plot-stats false'
'--levelling-para-compare-to-existing-model <logical>   e.g. --levelling-para-compare-to-existing-model true'
'--levelling-para-existing-model <path_to_file>         e.g. --levelling-para-existing-model /g/data/dg9/nd2979/Data/EXISTING_GEOID_MODELS/AGQG20221120.mat'
'--levelling-para-max-diff <value>                      e.g. --levelling-para-max-diff 0.15'
'--output-para-grids-name <path_to_folder>              e.g. --output-para-grids-name /g/data/dg9/nd2979/outputs/GridsNENSWgg2degTile/'
'--output-para-tiles-dir-name <path_to_folder>          e.g. --output-para-tiles-dir-name /g/data/dh8/outputs/ResidualTilesNENSWgg2degTile/'
'--output-para-plot-grids <logical>                     e.g. --output-para-plot-grids false'
'--output-para-plotsfolder <path_to_folder>             e.g. --output-para-plotsfolder /g/data/dh8/outputs/plots/22-Nov-2024NENSWgg2degTile'
'--keepawake <logical>                                  e.g. --keepawake true'
' '
'Geoscience Australia. Neda Darbeheshti on 24/11/2024'};

for i=1:length(str)
    disp(str{i})
end
return
