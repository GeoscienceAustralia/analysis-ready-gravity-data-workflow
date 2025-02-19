function geomGravDiff = mosaicTiles(GRID_PARA,DEM_PARA,OUTPUT_PARA,Lev,LongDEM,LatDEM,REFERENCE_GEOID_Zetai,GGM_Gi,GGM_Zetai,Coastline)
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
Grid_res_geoid=zeros(DEM_PARA.num_rows,DEM_PARA.num_cols);
Grid_res_geoid_err=zeros(DEM_PARA.num_rows,DEM_PARA.num_cols);
Grid_res_grav=zeros(DEM_PARA.num_rows,DEM_PARA.num_cols);
Grid_res_grav_Bouguer=zeros(DEM_PARA.num_rows,DEM_PARA.num_cols);
Grid_res_grav_err=zeros(DEM_PARA.num_rows,DEM_PARA.num_cols);
Weights=zeros(DEM_PARA.num_rows,DEM_PARA.num_cols);

ZDeg=mean(mean(REFERENCE_GEOID_Zetai(LongDEM,LatDEM)-GGM_Zetai(LongDEM,-LatDEM,LatDEM*0)));

resAGQG=REFERENCE_GEOID_Zetai(LongDEM,LatDEM)-GGM_Zetai(LongDEM,-LatDEM,LatDEM*0);

%for k=1:2
 for k=1:length(Files)

    Tile_Data=importdata([OUTPUT_PARA.Tiles_dir_name,Files(k).name]);
    Wf=Tile_Data.weights;
    Grid_res_geoid=Grid_res_geoid+(Wf).*reshape(Tile_Data.res_geoid,DEM_PARA.num_rows,DEM_PARA.num_cols);
    Grid_res_geoid_err=Grid_res_geoid_err+(Wf).*reshape(Tile_Data.pot_error,DEM_PARA.num_rows,DEM_PARA.num_cols);  
    Grid_res_grav=Grid_res_grav+(Wf).*reshape(Tile_Data.res_grav,DEM_PARA.num_rows,DEM_PARA.num_cols);  
    Grid_res_grav_Bouguer=Grid_res_grav_Bouguer+(Wf).*reshape(Tile_Data.res_grav_Bouguer,DEM_PARA.num_rows,DEM_PARA.num_cols);
    Grid_res_grav_err=Grid_res_grav_err+(Wf).*reshape(Tile_Data.grav_error,DEM_PARA.num_rows,DEM_PARA.num_cols);
    Weights=(Weights+(Wf));

    Grid_res_geoid_w=Grid_res_geoid./Weights;
    Grid_res_geoid_err_w=Grid_res_geoid_err./Weights;
    Grid_res_grav_w=Grid_res_grav./Weights;
    Grid_res_grav_Bouguer_w=Grid_res_grav_Bouguer./Weights;
    Grid_res_grav_err_w=Grid_res_grav_err./Weights;
    
    Geoid_temp=double(Grid_res_geoid_w+GGM_Zetai(LongDEM,-LatDEM,LatDEM*0));
       
    geoidLSCgriddedInterpolant=griddedInterpolant(LongDEM(end:-1:1,:)',LatDEM(end:-1:1,:)',Geoid_temp(end:-1:1,:)');
    
    geomGravDiff=Lev(:,3)-geoidLSCgriddedInterpolant(Lev(:,1),Lev(:,2));  
    
    AGQG_Vals_Lev=Lev(:,3)-REFERENCE_GEOID_Zetai(Lev(:,1),Lev(:,2)); 

end

if OUTPUT_PARA.PLOT_GRIDS
        plotMosaicTiles(Coastline,GRID_PARA,LongDEM,LatDEM,Grid_res_geoid_w,resAGQG,ZDeg,Lev,geomGravDiff, AGQG_Vals_Lev, ...
    Grid_res_geoid_err_w,Grid_res_grav_w,Grid_res_grav_Bouguer_w,Grid_res_grav_err_w,OUTPUT_PARA.plotsFolder)
end 

DisplayAreaStatistics(Coastline,GRID_PARA,LongDEM,LatDEM,Grid_res_geoid_w, ...
    Grid_res_geoid_err_w,OUTPUT_PARA)

% disp('Save mat files')
% 
% save([OUTPUT_PARA.Grids_name,'geomGravDiff',date,'.mat'],'geomGravDiff')
% 
% save([OUTPUT_PARA.Grids_name,'Grid_res_geoid_w',date,'.mat'],'Grid_res_geoid_w')
% 
% save([OUTPUT_PARA.Grids_name,'Grid_res_geoid_err_w',date,'.mat'],'Grid_res_geoid_err_w')
% 
% save([OUTPUT_PARA.Grids_name,'Grid_res_grav_w',date,'.mat'],'Grid_res_grav_w')
% 
% save([OUTPUT_PARA.Grids_name,'Grid_res_grav_err_w',date,'.mat'],'Grid_res_grav_err_w')
% 
% save([OUTPUT_PARA.Grids_name,'Grid_res_grav_Bouguer_w',date,'.mat'],'Grid_res_grav_Bouguer_w')

disp('Preparing final grids')

Grid_res_geoid_w(isnan(Grid_res_geoid_w))=0;
Grid_res_geoid_err_w(isnan(Grid_res_geoid_err_w))=-99999;
Grid_res_grav_w(isnan(Grid_res_grav_w))=0;
Grid_res_grav_err_w(isnan(Grid_res_grav_err_w))=-99999;
Grid_res_grav_Bouguer_w(isnan(Grid_res_grav_Bouguer_w))=NaN;

final_quasigeoid_model=double(Grid_res_geoid_w+ZDeg+GGM_Zetai(LongDEM,-LatDEM,LatDEM*0));
final_quasigeoid_model_err=Grid_res_geoid_err_w;

final_freeair_gravity_model=double(Grid_res_grav_w+GGM_Gi(LongDEM,-LatDEM,LatDEM*0));
final_freeair_gravity_model_err=Grid_res_grav_err_w;

final_Bouguer_gravity_model=double(Grid_res_grav_Bouguer_w+GGM_Gi(LongDEM,-LatDEM,LatDEM*0));
final_Bouguer_gravity_model(isnan(Grid_res_grav_Bouguer_w))=-99999;

[Long_out,Lat_out]=meshgrid(93:1/60:174-1/60,-8:-1/60:-61+1/60);

% resample on the existing grid.

resamplegeoid=griddata(LongDEM,LatDEM,final_quasigeoid_model,Long_out,Lat_out);
resamplegeoid_err=griddata(LongDEM,LatDEM,final_quasigeoid_model_err,Long_out,Lat_out);
resamplegravity=griddata(LongDEM,LatDEM,final_freeair_gravity_model,Long_out,Lat_out);
resamplegravity_bouguer=griddata(LongDEM,LatDEM,final_Bouguer_gravity_model,Long_out,Lat_out);
resamplegravity_err=griddata(LongDEM,LatDEM,final_freeair_gravity_model_err,Long_out,Lat_out);

disp('Saving to geotiffs')

% GCSE_GRS1980 =	4019 from http://geotiff.maptools.org/spec/geotiff6.html#6.2

% Define the bounding box
%bbox = [93, -61 - 1/60; 174 - 1/60, -8];
bbox = [93, -61 ; 174 , -8];
% Set the GeoTIFF options
option.GeographicTypeGeoKey = 4019;
option.GTModelTypeGeoKey  = 2;
option.GTRasterTypeGeoKey = 1;

% Specify bit depth
bit_depth = 32;
formattedDate = datestr(date, 'yyyymmdd');

datwrite(Long_out,Lat_out,resamplegeoid,OUTPUT_PARA,formattedDate);

geotiffwrite([OUTPUT_PARA.Grids_name,'AGQG_',formattedDate,'.tif'],bbox, resamplegeoid, bit_depth, option);
geotiffwrite([OUTPUT_PARA.Grids_name,'AGQG_1SigUncert_',formattedDate,'.tif'],bbox, resamplegeoid_err, bit_depth, option);
geotiffwrite([OUTPUT_PARA.Grids_name,'AGQG_Free_Air_Anomaly_',formattedDate,'.tif'],bbox, resamplegravity, bit_depth, option);
geotiffwrite([OUTPUT_PARA.Grids_name,'AGQG_Bouguer_Anomaly_',formattedDate,'.tif'],bbox, resamplegravity_bouguer, bit_depth, option);
geotiffwrite([OUTPUT_PARA.Grids_name,'AGQG_Gravity_1SigUncert_',formattedDate,'.tif'],bbox, resamplegravity_err, bit_depth, option);

end






    










