
GRID_PARA.STEP=0.5;
GRID_PARA.MINLONG=144.5;%93;
GRID_PARA.MAXLONG=150;%173;
GRID_PARA.MINLAT=-39;%-59;
GRID_PARA.MAXLAT=-36.5;%-9;

lat_range = GRID_PARA.MAXLAT:-GRID_PARA.STEP:GRID_PARA.MINLAT;
long_range = GRID_PARA.MINLONG:GRID_PARA.STEP:GRID_PARA.MAXLONG;

n_lat = length(lat_range);
n_long = length(long_range);

% Total number of tiles, Justy' formula
n_total = (n_lat-1) * (n_long-1);