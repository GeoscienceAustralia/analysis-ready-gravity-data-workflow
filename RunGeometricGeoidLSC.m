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
% NSW=[139.5 154 -38 -27.5];
GRID_PARA.MINLONG=140;
GRID_PARA.MAXLONG=154;
GRID_PARA.MINLAT=-38;
GRID_PARA.MAXLAT=-29;
%% DEM data - N.B. the dem is used to specify the grid nodes.
DEM_PARA.filename='Data\DEM\AUSDEM1min.xyz';
DEM_PARA.num_cols=4861;
DEM_PARA.num_rows=3181;
%% Gravity data
% First run ./Data/GRAVITY/XXXX/PrepareGravity_XXXXX.m
% And then /Data/GRAVITY/Combine_Gravity_Data.m
% this collates all of the gravity and position data into one matlab array.
GRAV_PARA.filename='Data/processedData/Terrestrial_Gravity.mat';%'Data\GRAVITY\Gravity.mat';
GRAV_PARA.filename1=[];%'Data\GRAVITY\Xcalibur_Gravity.mat';
GRAV_PARA.TypeB=1;% This is a Type B uncertainty value (in mGal) which is added to the uncertainty values.
GRAV_PARA.Grav_Faye_TypeB=3;
%% Gravity Gradiometry data
% Add notes here
GRAV_GRAD_PARA.filename='Data\GRAVITY_GRAD\Xcalibur_FVD_GDD.mat';
GRAV_GRAD_PARA.TypeB=10^(-5);% This is a Type B uncertainty value (in mGal/m) which is added to the uncertainty values.
GRAV_GRAD_PARA.avail=false;
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
COV_PARA.COVPlot=false;% true plots progress, false turns this off.
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
LEVELLING_PARA.filename='Data\GPS_LEVELLING\Lev_NSW_NG.mat';%'Data/GPS_LEVELLING/Lev_CARS.mat';% The format of these data needs to be an array with rows [Long,Lat,h-H].
LEVELLING_PARA.Plot_Stats=true;% If true, the levelling data are compared to the geoid as its computed.
LEVELLING_PARA.Compare_To_Existing_Model=true;% If true, the levelling data are also compared to another existing geoid as its computed.
LEVELLING_PARA.Existing_Model='Data\EXISTING_GEOID_MODELS\AGQG20221120.mat';% File location of the existing model.
LEVELLING_PARA.max_diff=0.15;% Threshold for an outlier with the GNSS-levelling
%% Output
OUTPUT_PARA.Grids_name='D:\GAbackup\Outputs\GridsVicNSWterr/';
OUTPUT_PARA.Tiles_dir_name='D:/GAbackup/Outputs/ResidualTilesvicNSWterr/';
OUTPUT_PARA.PLOT_GRIDS=false;% A gridded solution is plotted and output as well as the tiles.
OUTPUT_PARA.plotsFolder=['D:/GAbackup/Outputs/plots/',date,'vicNSWterr'];
% Keep the computer awake
keepawake=true;% Setting this to true wiggles the mouse every so often so the compute doesnt go to sleep.
%% Run the LSC code

disp('1/4 ..........................importAndFormatData is running ')
[Gravo,gravGradFiltered,DEM_data,ZDEM_griddedInterpolant,LongDEM,LatDEM,...
 GGM_Gravity_griddedInterpolant,GGM_Zeta_griddedInterpolant,Lev,...
 REFERENCE_Zeta_griddedInterpolant,GRID_REF,Coastline]=importAndFormatData...
 (GRID_PARA,DEM_PARA,GRAV_PARA,Topo_PARA,COAST_PARA,LEVELLING_PARA,GGM_PARA,GRAV_GRAD_PARA);

% disp('4/4 ..........................mosaicTiles is running')
% Vals_Lev = MosaicTiles(GRID_PARA,DEM_PARA,OUTPUT_PARA,Lev,LongDEM,LatDEM, ...
%     REFERENCE_Zeta_griddedInterpolant,GGM_Gravity_griddedInterpolant,GGM_Zeta_griddedInterpolant,Coastline);

%save([OUTPUT_PARA.Grids_name,'levellingLSCterrGeoids',date,'.mat'],'Vals_Lev')
Vals_Lev=importdata([OUTPUT_PARA.Grids_name,'levellingLSCterrGeoids18-Jun-2024.mat']);

%% Remove a tiled plane so the signal is zero mean for the LSC
coeffs=[Lev(:,1)-mean(Lev(:,1)),Lev(:,2)-mean(Lev(:,2)),ones(size(Lev(:,2)))]\Vals_Lev;
Vals_Lev_detrended=Vals_Lev-[Lev(:,1)-mean(Lev(:,1)),Lev(:,2)-mean(Lev(:,2)),ones(size(Lev(:,2)))]*coeffs;

% plot up the residuals
figure(8)
clf
hold on
scatter(Lev(:,1),Lev(:,2),15,Vals_Lev_detrended,'filled')
% scatter(Lev(:,1),Lev(:,2),15,Vals_Lev,'o')
colorbar
for jj=1:9
    plot(Coastline.long{jj},Coastline.lat{jj},'k')
end
colorbar
caxis([-0.25 0.25])
colormap(jet)
ax=gca;
ax.PlotBoxAspectRatio =[1 abs(cos(mean(LatDEM(:,1))*pi/180)) 1];
axis tight 
axis([148.5 150.5 -34.2 -32.5])

%% Compute auto covariance
% Estimate empirical covariance
Out=computeSphericalEmpiricalCovariance(Lev(:,1),Lev(:,2),Vals_Lev_detrended,1);
disp('Pre-computing covariance function')
[EMP_COV,psi,sigma2,l]=fitEmpiricalCovarianceB(Out(:,1),Out(:,2));
%% Compute the LSC linear functionals
for k=1:length(Lev(:,1))
sinPSIon2=sqrt(sin((Lev(k,2)*pi/180-Lev(:,2)*pi/180)/2).^2+(sin((Lev(k,1)*pi/180-Lev(:,1)*pi/180)/2).^2).*cos(Lev(k,2)*pi/180).*cos(Lev(:,2)*pi/180));
psi=2*asin(sinPSIon2');
ACOV_tt(k,:)=sigma2*exp(-(psi.^2)/(2*l.^2));
end

% inverse of auto covariance matrix
 inverseCovarianceMatrix=(ACOV_tt+0.000025*eye(size(ACOV_tt)))\eye(size(ACOV_tt));


 sol=inverseCovarianceMatrix*(Vals_Lev_detrended);



%sol=inv(ACOV_tt+0.000025*eye(size(ACOV_tt)))*(Vals_Lev_detrended);

%sol=(ACOV_tt+0.000025*eye(size(ACOV_tt)))\(Vals_Lev_detrended);
%%  Now do the interpolation one row of latitude at a time.
LSC_sol=LongDEM*0;
LSC_solrt=LSC_sol;
for ki=1:length(LongDEM(:,1))
ACOV_tt_dem=[];
for k=1:length(Lev(:,1))
sinPSIon2=sqrt(sin((Lev(k,2)*pi/180-LatDEM(ki,:)*pi/180)/2).^2+(sin((Lev(k,1)*pi/180-LongDEM(ki,:)*pi/180)/2).^2).*cos(Lev(k,2)*pi/180).*cos(LatDEM(ki,:)*pi/180));
psi=2*asin(sinPSIon2);
ACOV_tt_dem(k,:)=sigma2*exp(-(psi.^2)/(2*l.^2));
end
ACOV_tt_dem=ACOV_tt_dem';

LSC_sol(ki,:)=ACOV_tt_dem*sol;
disp(ki)
if mod((length(LongDEM(1,:))-ki),10)==0
    figure(9)
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),LSC_sol)
    scatter(Lev(:,1),Lev(:,2),15,Vals_Lev_detrended,'filled')
    title(length(LongDEM(1,:))-ki)
    for jj=1:9
    plot(Coastline.long{jj},Coastline.lat{jj},'k')
    end
    colorbar
    colormap(jet)
    axis([148.5 150.5 -34.2 -32.5])   
    drawnow

LSC_solrt(:)=LSC_sol(:)+[LongDEM(:)-mean(Lev(:,1)),LatDEM(:)-mean(Lev(:,2)),ones(size(LongDEM(:)))]*coeffs;

    figure(10)
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),LSC_solrt)
    scatter(Lev(:,1),Lev(:,2),15,Vals_Lev,'filled')
%     plot(Lev(:,1),Lev(:,2),'ko')
    title(length(LongDEM(1,:))-ki)
    for jj=1:9
    plot(Coastline.long{jj},Coastline.lat{jj},'k')
    end
    colorbar
    caxis([min(Vals_Lev) max(Vals_Lev)])
    colormap(jet)
    axis([148.5 150.5 -34.2 -32.5])   
    drawnow

end
end
%% Add code to restore the tilt, then add back the reference geoid model.
%.
%..
%...

%% functions for the LSC



