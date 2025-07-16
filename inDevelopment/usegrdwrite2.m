filename = 'Data/EXISTING_GEOID_MODELS/AUSGeoid2020_20180201_win.dat'; % change to your file name/path

fid = fopen(filename, 'r');
fgetl(fid); % skip header

line = fgetl(fid);
disp(['Sample line: ', line]);
% Relax spacing by using spaces, and try:
data = textscan(fid, 'GEO %f S %d %d %f E %d %d %f %f %f');
fclose(fid);

lat_deg  = double(data{2});  % degrees
lat_min  = double(data{3});  % minutes
% Convert latitude to decimal degrees (South → negative)
latitude_decimal = -(lat_deg + lat_min / 60);

lon_deg  = double(data{5});  % degrees
lon_min  = double(data{6});  % minutes
% Convert latitude to decimal degrees (South → negative)
longitude_decimal = (lon_deg + lon_min / 60);

grdwrite2(longitude_decimal,latitude_decimal,double(data{1}),'Data/EXISTING_GEOID_MODELS/AUSGeoid2020_20180201_win.grd')


