GRID_PARA.STEP=0.5;
GRID_PARA.MINLONG=94;
GRID_PARA.MAXLONG=173;
GRID_PARA.MINLAT=-60;
GRID_PARA.MAXLAT=-9;

lat_range = GRID_PARA.MAXLAT:-GRID_PARA.STEP:GRID_PARA.MINLAT;
long_range = GRID_PARA.MINLONG:GRID_PARA.STEP:GRID_PARA.MAXLONG;

n_lat = length(lat_range);
n_long = length(long_range);

% Total number of tiles
n_total = n_lat * n_long;