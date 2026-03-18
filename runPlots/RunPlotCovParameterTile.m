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
% Add the path to the function files.
addpath('functions');
% Grid/Tiling Parameters
GRID_PARA.buffer=1;% degs. The x/y extent to extract data around the tile centre.
GRID_PARA.MINLONG=150;
GRID_PARA.MAXLONG=151;
GRID_PARA.MINLAT=-28;
GRID_PARA.MAXLAT=-27;

COAST_PARA.filename='Data\COASTLINE\CoastAus.mat';
% Import coastline
Coastline=importdata(COAST_PARA.filename);

OUTPUT_PARA.Tiles_dir_name='Outputs/ResidualTilesNQueen150longPar/';
Files=dir(OUTPUT_PARA.Tiles_dir_name);
Files(1:2)=[];

for k=1:length(Files)
    Tile_Data=importdata([OUTPUT_PARA.Tiles_dir_name,Files(k).name]);
    % Access the name of the file from the struct
    filename = Files(k).name;
    % Find the index of the first underscore
    underscore_index = strfind(filename, '_');
    m_index = strfind(filename, 'm');
    % Extract characters from the fourth position to the underscore
    filenameLong(k) = str2double(filename(5:underscore_index - 1));
    filenameLat(k) =str2double(filename(underscore_index+1:m_index-2));

    COV_PARA_RTM_A(k)=Tile_Data.COV_PARA_RTM.A;
    COV_PARA_RTM_B(k)=Tile_Data.COV_PARA_RTM.B;
    COV_PARA_RTM_M(k)=Tile_Data.COV_PARA_RTM.M;
    COV_PARA_Faye_A(k)=Tile_Data.COV_PARA_Faye.A;
    COV_PARA_Faye_B(k)=Tile_Data.COV_PARA_Faye.B;
    COV_PARA_Faye_M(k)=Tile_Data.COV_PARA_Faye.M;
end




% common variables for plotting
axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;

% plots
figure('Name','MosaicTiles','NumberTitle','off'); 
hold on
scatter (filenameLong,-filenameLat,55,COV_PARA_RTM_A,'filled')
customizeMap('COV PARA RTM A','',Coastline,axisLimits)
%caxis([0 220])

figure('Name','MosaicTiles','NumberTitle','off'); 
hold on
scatter (filenameLong,-filenameLat,55,COV_PARA_RTM_B,'filled')
customizeMap('COV PARA RTM B','',Coastline,axisLimits)

figure('Name','MosaicTiles','NumberTitle','off'); 
hold on
scatter (filenameLong,-filenameLat,55,COV_PARA_RTM_M,'filled')
customizeMap('COV PARA RTM M','',Coastline,axisLimits)
%caxis([0 220])


save(['covParameters',date,'.mat'],'COV_PARA_RTM_A','COV_PARA_RTM_B','COV_PARA_RTM_M',...,
                                   'COV_PARA_Faye_A','COV_PARA_Faye_B','COV_PARA_Faye_M',...,
                                   'filenameLong','filenameLat')


 

