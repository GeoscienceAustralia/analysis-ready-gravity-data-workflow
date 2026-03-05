function geomGravDiff = mosaicTilesParallel( ...
    GRID_PARA, DEMpara, OUTPUT_PARA, Lev, ...
    LongDEM, LatDEM, REFERENCE_GEOID_Zetai, ...
    GGM_Gi, GGM_Zetai, Coastline)

disp('Tiling...')

Files = dir(OUTPUT_PARA.Tiles_dir_name);
Files(1:2) = [];
nFiles = numel(Files);

nR = DEMpara.num_rows;
nC = DEMpara.num_cols;

% ---- Reduction accumulators ----
sum_geoid      = zeros(nR,nC);
sum_geoid_err  = zeros(nR,nC);
sum_grav       = zeros(nR,nC);
sum_grav_boug  = zeros(nR,nC);
sum_grav_err   = zeros(nR,nC);
sum_weights    = zeros(nR,nC);

tilesDir = OUTPUT_PARA.Tiles_dir_name;

% ---- Parallel accumulation ----
parfor k = 1:nFiles

    Tile_Data = importdata(fullfile(tilesDir, Files(k).name));
    Wf = Tile_Data.weights;

    sum_geoid      = sum_geoid      + Wf .* reshape(Tile_Data.res_geoid,        nR, nC);
    sum_geoid_err  = sum_geoid_err  + Wf .* reshape(Tile_Data.pot_error,        nR, nC);
    sum_grav       = sum_grav       + Wf .* reshape(Tile_Data.res_grav,         nR, nC);
    sum_grav_boug  = sum_grav_boug  + Wf .* reshape(Tile_Data.res_grav_Bouguer, nR, nC);
    sum_grav_err   = sum_grav_err   + Wf .* reshape(Tile_Data.grav_error,       nR, nC);
    sum_weights    = sum_weights    + Wf;

end

% ---- Final weighted grids ----
Weights = sum_weights;
Weights(Weights == 0) = NaN;

Grid_res_geoid_w        = sum_geoid      ./ Weights;
Grid_res_geoid_err_w    = sum_geoid_err  ./ Weights;
Grid_res_grav_w         = sum_grav       ./ Weights;
Grid_res_grav_Bouguer_w = sum_grav_boug  ./ Weights;
Grid_res_grav_err_w     = sum_grav_err   ./ Weights;

% ---- Reference offsets ----
ZDeg = mean(mean( ...
    REFERENCE_GEOID_Zetai(LongDEM,LatDEM) - ...
    GGM_Zetai(LongDEM,-LatDEM,LatDEM*0)));

resAGQG = REFERENCE_GEOID_Zetai(LongDEM,LatDEM) - ...
          GGM_Zetai(LongDEM,-LatDEM,LatDEM*0);

% ---- Interpolant (ONCE) ----
Geoid_temp = double(Grid_res_geoid_w + ...
    GGM_Zetai(LongDEM,-LatDEM,LatDEM*0));

geoidLSCgriddedInterpolant = griddedInterpolant( ...
    LongDEM(end:-1:1,:)', ...
    LatDEM(end:-1:1,:)', ...
    Geoid_temp(end:-1:1,:)' );

geomGravDiff = Lev(:,3) - ...
    geoidLSCgriddedInterpolant(Lev(:,1), Lev(:,2));

AGQG_Vals_Lev = Lev(:,3) - ...
    REFERENCE_GEOID_Zetai(Lev(:,1), Lev(:,2));

% ---- Plotting ----
if OUTPUT_PARA.PLOT_GRIDS
    plotMosaicTiles( ...
        Coastline, GRID_PARA, LongDEM, LatDEM, ...
        Grid_res_geoid_w, resAGQG, ZDeg, Lev, ...
        geomGravDiff, AGQG_Vals_Lev, ...
        Grid_res_geoid_err_w, Grid_res_grav_w, ...
        Grid_res_grav_Bouguer_w, Grid_res_grav_err_w, ...
        OUTPUT_PARA.plotsFolder)
end

% ---- Statistics ----
DisplayAreaStatistics( ...
    Coastline, GRID_PARA, LongDEM, LatDEM, ...
    Grid_res_geoid_w, Grid_res_geoid_err_w, OUTPUT_PARA)

disp('Save mat files')

save([OUTPUT_PARA.Grids_name,'geomGravDiff',date,'.mat'],'geomGravDiff')
save([OUTPUT_PARA.Grids_name,'Grid_res_geoid_w',date,'.mat'],'Grid_res_geoid_w')
save([OUTPUT_PARA.Grids_name,'Grid_res_geoid_err_w',date,'.mat'],'Grid_res_geoid_err_w')
save([OUTPUT_PARA.Grids_name,'Grid_res_grav_w',date,'.mat'],'Grid_res_grav_w')
save([OUTPUT_PARA.Grids_name,'Grid_res_grav_err_w',date,'.mat'],'Grid_res_grav_err_w')
save([OUTPUT_PARA.Grids_name,'Grid_res_grav_Bouguer_w',date,'.mat'],'Grid_res_grav_Bouguer_w')

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
