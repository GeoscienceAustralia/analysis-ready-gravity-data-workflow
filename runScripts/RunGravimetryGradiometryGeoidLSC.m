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
%
close all
clear 
% Turn off the irritating warnings.
warning off
% Add the path to the function files.
addpath('functions');
% run myBlanketFunction to define all parameters 
load('allParameters.mat')

disp('1/4 ..........................importAndFormatData is running ')
[Gravo,gravGradFiltered,DEM_data,ZDEM_griddedInterpolant,LongDEM,LatDEM,...
 GGM_Gravity_griddedInterpolant,GGM_Zeta_griddedInterpolant,Lev,...
 REFERENCE_Zeta_griddedInterpolant,GRID_REF,Coastline,DEM_PARA]=importAndFormatData...
 (GRID_PARA,DEM_PARA,GRAV_PARA,Topo_PARA,COAST_PARA,LEVELLING_PARA,GGM_PARA,GRAV_GRAD_PARA);

if OUTPUT_PARA.PLOT_GRIDS
     plotInputData(Gravo,Coastline,GRID_PARA)
end 

if GRAV_PARA.altimetry_weighting 
     weightAltimetry(Gravo,Coastline,GRID_PARA,OUTPUT_PARA)
end

disp('2/4 ..........................computeTerrainEffect is running')
[fullTopoCorrectedGravityPoint,longwaveTopo_griddedInterpolant,fullTopo_griddedInterpolant,fullTopoCorrectedGravityGradient]=computeFullTerrainEffects(GRID_PARA, ...
    Topo_PARA,Gravo,gravGradFiltered,GGM_Gravity_griddedInterpolant,DEM_data,ZDEM_griddedInterpolant, ...
    LongDEM,LatDEM,Coastline,OUTPUT_PARA.plotsFolder);

save([OUTPUT_PARA.Grids_name,'TerrainEffects',date,'.mat'],'fullTopoCorrectedGravityPoint','longwaveTopo_griddedInterpolant','fullTopo_griddedInterpolant','fullTopoCorrectedGravityGradient')

disp('3/4 ..........................computeGravimetryGradiometryLSC is running')
computeGravimetryGradiometryLSC(GRID_PARA,COV_PARA,DEM_PARA,GRAV_PARA,GRAV_GRAD_PARA,OUTPUT_PARA,GRID_REF,fullTopoCorrectedGravityPoint,fullTopoCorrectedGravityGradient, ...
    GGM_Gravity_griddedInterpolant,ZDEM_griddedInterpolant,fullTopo_griddedInterpolant, ...
    longwaveTopo_griddedInterpolant,Topo_PARA.Density)

% uncomment if there is a mat file for TerrainEffects
% TE = importdata([OUTPUT_PARA.Grids_name,'TerrainEffects29-Oct-2024.mat']);
% disp('3/4 ..........................computeGravimetryGradiometryLSC is running')
% computeGravimetryGradiometryLSC(GRID_PARA,COV_PARA,DEM_PARA,GRAV_PARA,GRAV_GRAD_PARA,OUTPUT_PARA,GRID_REF,TE.fullTopoCorrectedGravityPoint,TE.fullTopoCorrectedGravityGradient, ...
%     GGM_Gravity_griddedInterpolant,ZDEM_griddedInterpolant,TE.fullTopo_griddedInterpolant, ...
%     TE.longwaveTopo_griddedInterpolant,Topo_PARA.Density)

disp('4/4 ..........................mosaicTiles is running')
geomGravGeoidDiff = mosaicTiles(GRID_PARA,DEM_PARA,OUTPUT_PARA,Lev,LongDEM,LatDEM, ...
    REFERENCE_Zeta_griddedInterpolant,GGM_Gravity_griddedInterpolant,GGM_Zeta_griddedInterpolant,Coastline);



