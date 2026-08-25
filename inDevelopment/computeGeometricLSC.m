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

<<<<<<< HEAD
 plotLevellingData(Lev,Lev(:,3),'h - H','m',OUTPUT_PARA.plotsFolder)
=======
plotLevellingData(Lev(:,1),Lev(:,2), Lev(:,3), 'h-H', OUTPUT_PARA.plotsFolder)
>>>>>>> 9878bc5efb7e47c1338044d5c71ab46c797435b9

% read final matfiles

dateCreated ='23-Mar-2026';

load([OUTPUT_PARA.Grids_name,'Grid_res_geoid_w',dateCreated,'.mat'])

ZDeg=mean(mean(REFERENCE_Zeta_griddedInterpolant(LongDEM,LatDEM)-GGM_Zeta_griddedInterpolant(LongDEM,-LatDEM,LatDEM*0)));

resAGQG=REFERENCE_Zeta_griddedInterpolant(LongDEM,LatDEM)-GGM_Zeta_griddedInterpolant(LongDEM,-LatDEM,LatDEM*0);

Geoid_temp=double(Grid_res_geoid_w+GGM_Zeta_griddedInterpolant(LongDEM,-LatDEM,LatDEM*0));
       
geoidLSCgriddedInterpolant=griddedInterpolant(LongDEM(end:-1:1,:)',LatDEM(end:-1:1,:)',Geoid_temp(end:-1:1,:)');
    
geomGravDiff=Lev(:,3)-geoidLSCgriddedInterpolant(Lev(:,1),Lev(:,2));  

%geomGravDiff2022=Lev(:,3)-REFERENCE_Zeta_griddedInterpolant(Lev(:,1),Lev(:,2));

<<<<<<< HEAD
plotLevellingData(Lev,geomGravDiff,'h - H - AGQG','m',OUTPUT_PARA.plotsFolder)

if LEVELLING_PARA.Plot_Stats
   plotGPSlevelling(Coastline,GRID_PARA,Lev,geomGravDiff,geomGravDiff2022,OUTPUT_PARA.plotsFolder)
end 
=======
plotLevellingData(Lev(:,1),Lev(:,2), geomGravDiff,[],'h-H-AGQG', OUTPUT_PARA.plotsFolder)
>>>>>>> 9878bc5efb7e47c1338044d5c71ab46c797435b9

% Remove a tiled plane so the signal is zero mean for the LSC

[plane_params, res_detrended] = remove_tilted_plane_geo(Lev(:,2), Lev(:,3), geomGravDiff); % from pyhton

plotLevellingData(Lev(:,1),Lev(:,2), res_detrended,[-0.3 0.3],'h-H-AGQG Detrended Python', OUTPUT_PARA.plotsFolder)

% Construct the matrix for linear trend removal
trendMatrix = [Lev(:,1) - mean(Lev(:,1)), Lev(:,2) - mean(Lev(:,2)), ones(size(Lev(:,2)))];

% Calculate the coefficients of the best-fit plane
trendCoefficients = trendMatrix \ geomGravDiff;
%trendCoefficients2022 = trendMatrix \ geomGravDiff2022;

% Remove the planar trend to obtain zero-mean data
geomGravGeoidDiffDetrended = geomGravDiff - trendMatrix * trendCoefficients;
%geomGravGeoidDiff2022Detrended = geomGravDiff2022 - trendMatrix * trendCoefficients2022;

plotLevellingData(Lev(:,1),Lev(:,2), geomGravGeoidDiffDetrended,[-0.3 0.3],'h-H-AGQG Detrended', OUTPUT_PARA.plotsFolder)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
disp('computing covariance functions')

covarianceInfo=computeSphericalEmpiricalCovariance(Lev(:,1),Lev(:,2),geomGravGeoidDiffDetrended,1);

[sigma2, correlationL, Cfit] = fitGaussianCovariance(covarianceInfo(:,1),covarianceInfo(:,2),OUTPUT_PARA.plotsFolder, ...
    'anchor_sigma2', false, ...
    'do_plot', true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert degrees to radians
longitudeLevRadian = deg2rad (Lev(:,1));
latitudeLevRadian = deg2rad (Lev(:,2));

% initialize covariance matrix
ACOVtt = zeros(length(Lev(:,1)),length(Lev(:,1)));
for lonCounter=1:length(Lev(:,1))
haversineDistance=haversine(latitudeLevRadian(lonCounter), longitudeLevRadian(lonCounter),latitudeLevRadian(:), longitudeLevRadian(:));
ACOVtt(lonCounter,:)=sigma2*exp(-(haversineDistance.^2)/(2*correlationL.^2));
end

% Plot covariance matrix
plotCovarianceMatrix(haversineDistance, ACOVtt,'m^2','Gaussian Covariance', OUTPUT_PARA.plotsFolder)

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psi_vals = linspace(0, pi/6, 10000);

C_tt_interp = build_covariance_interpolator_mixed_gaussian(psi_vals, sigma2, correlationL);

C_tt = evaluate_covariance_block_mixed_gaussian(Lev(:,2), Lev(:,1), Lev(:,2) , Lev(:,1), C_tt_interp);

% Plot covariance matrix
plotCovarianceMatrix(haversineDistance, C_tt,'m^2','Gaussian CovariancePyhton', OUTPUT_PARA.plotsFolder)

% Solve linear system:
% rhs = (C_tt + 0.0025 * I) \ res_detrended
n = numel(geomGravGeoidDiffDetrended);
rhs = (C_tt + 0.0025 * eye(n)) \ geomGravGeoidDiffDetrended;

% Print standard deviations
fprintf('std(res_detrended - C_tt * rhs) = %.6f\n', ...
        std(geomGravGeoidDiffDetrended - C_tt * rhs));

fprintf('std(res_detrended) = %.6f\n', ...
        std(geomGravGeoidDiffDetrended));

% Plot residuals after LSC
plotLevellingData(Lev(:,1),Lev(:,2), geomGravGeoidDiffDetrended - C_tt * rhs,[-0.3 0.3],'LSC Residual', OUTPUT_PARA.plotsFolder)

% Initialise output grid (same shape as lat_grid_potential_restore)
tt_gp = LatDEM * 0;

% Grid size
[n_rows, n_cols] = size(LongDEM);

% Choose a block size (same meaning as Python)
block_size = 1;   % adjust as needed

for i0 = 1:block_size:n_rows
    fprintf('%d\n', i0);   % equivalent to print(i0)

    i1 = min(i0 + block_size - 1, n_rows);

    % Take a block of rows and flatten to 1D for the covariance call
    lat_block = reshape(LatDEM(i0:i1, :), [], 1);
    lon_block = reshape(LongDEM(i0:i1, :), [], 1);

    % Evaluate covariance block
    Ctt_block = evaluate_covariance_block_mixed_gaussian( ...
                    lat_block, lon_block, ...
                    Lev(:,2), Lev(:,1), ...
                    C_tt_interp );

    % Multiply by solution vector
    tt_block_flat = Ctt_block * rhs;

    % Reshape back to (n_block_rows x n_cols) and store
    tt_gp(i0:i1, :) = reshape(tt_block_flat, i1 - i0 + 1, n_cols);
end

[tt_gp_restored, plane_grid] = restore_tilted_plane_geo( ...
    LongDEM, ...
    LatDEM, ...
    tt_gp, ...
    plane_params);

plotGridedGeometriCorrection( tt_gp, LongDEM, LatDEM,[-0.3 0.3], ...
    GRID_PARA, OUTPUT_PARA, Coastline,'Geometric Correction using LSC')

plotGridedGeometriCorrection( tt_gp_restored, LongDEM, LatDEM,[], ...
    GRID_PARA, OUTPUT_PARA, Coastline,'Geometric Correction with trend using LSC')

% figure
% imagesc(Ctt_block)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use matlab generic interpolator
 % Extract coordinates
lon = Lev(:,1);
lat = Lev(:,2);
z   = geomGravGeoidDiffDetrended;

% --------------------------------------------------------------
% Remove invalid values
idx = isfinite(lon) & isfinite(lat) & isfinite(z);
lon = lon(idx);
lat = lat(idx);
z   = z(idx);

% --------------------------------------------------------------
% Create interpolant
F = scatteredInterpolant( ...
    lon, lat, z, ...
    'natural', ...   % interpolation
    'none');         % no extrapolation

% Evaluate on DEM grid
Zq = F(LongDEM, LatDEM);

plotGridedGeometriCorrection( Zq, LongDEM, LatDEM,[-0.3 0.3], ...
    GRID_PARA, OUTPUT_PARA, Coastline,'Geometric Correction using MATLAB')





