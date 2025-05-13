function Grid_res_grav_w = mosaicTiles4MolodenskyG1(GRID_PARA,DEMpara,OUTPUT_PARA,LongDEM,LatDEM,Coastline)
% Run this function to produce the final geoid model. It will:
% - Collate all tiles in the OUTPUT_PARA.Tiles_dir_name directory.
% - Add back the GGM.
% - Output gravity, geoid, and respective error values as TIFF files.
% The geoid model is aligned with the same zero degree term as for AGQG2017.
% Input:  GRID_PARA = grid/tiling parameters such as extent(MINLONG,MAXLONG,MINLAT,MAXLAT), buffer, STEP    
%         DEM_PARA=   DEM data such as filename,num_cols,num_rows
%         OUTPUT_PARA= output parameters such as Grids_name,Tiles_dir_name
%         LongDEM= matrix of longitudes for DEM from Import_And_Format_Data_Sets_Reg
%         LatDEM=  matrix of latitudes for DEM from Import_And_Format_Data_Sets_Reg
%         GGM_Gravity_griddedInterpolant = GGM gravity griddedInterpolant 
%         GGM_Zeta_griddedInterpolant = GGM height anomaly griddedInterpolant 
%         ZDEM_griddedInterpolant = DEM elevation griddedInterpolant
%         Lev= from Import_And_Format_Data_Sets_Reg?
%         REFERENCE_GEOID_Zetai=1*1 griddedInterpolant from Import_And_Format_Data_Sets_Reg?
%         Coastline= 1*1 struct from Import_And_Format_Data_Sets_Reg
%
% Output: geomGravDiff = a vector of differences between geometric and gravimetric geoid at GPS leveling points
%
% Example: see RunGeoidLSC
%
% Main functions
% - 
% Other functions
% - 
% - 
%
% Written by Jack McCubbine
% Last updated by Neda Darbeheshti
% Geoscience Australia, 2024-03.

% Import tiles and combine
disp('Tiling...')
Files=dir(OUTPUT_PARA.Tiles_dir_name);
Files(1:2)=[];

Grid_res_grav=zeros(DEMpara.num_rows,DEMpara.num_cols);
Weights=zeros(DEMpara.num_rows,DEMpara.num_cols);
%for k=1:2
 for k=1:length(Files)

    Tile_Data=importdata([OUTPUT_PARA.Tiles_dir_name,Files(k).name]);
    Wf=Tile_Data.weights;
    
    Grid_res_grav=Grid_res_grav+(Wf).*reshape(Tile_Data.res_grav,DEMpara.num_rows,DEMpara.num_cols);  
    Weights=(Weights+(Wf));
    Grid_res_grav_w=Grid_res_grav./Weights;
 end

 if OUTPUT_PARA.PLOT_GRIDS
    % common variables for plotting
    axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
    axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
    axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
    axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
    axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;
    % plot residualGravityWeighted
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_grav_w)
    customizeMap('Residual Free Air Gravity Weighted','mGal',Coastline,axisLimits)
    saveas(gcf,[OUTPUT_PARA.plotsFolder,'MosaicTiles','residualFreeAirGravityWeighted','.png'])
 end 

disp('Save mat files')

save([OUTPUT_PARA.Grids_name,'Grid_res_grav_w',date,'.mat'],'Grid_res_grav_w')

%save([OUTPUT_PARA.Grids_name,'Grid_res_grav_err_w',date,'.mat'],'Grid_res_grav_err_w')

end






    










