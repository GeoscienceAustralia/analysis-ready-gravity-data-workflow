% Path to your input dat file
inputFile ='Data/EXISTING_GEOID_MODELS/AUSGeoid2020_20180201_win.dat';

% Open the file for reading
fid = fopen(inputFile, 'r');

% Read lines (skip header)
data = textscan(fid, 'GEO %f S%d %f %f E%d %f %f %f %f', 'HeaderLines', 1);
fclose(fid);

% Extract and cast all to double
geoid      = double(data{1});
lat_deg    = double(data{2});
lat_min    = double(data{3});
lat_sec    = double(data{4});
lon_deg    = double(data{5});
lon_min    = double(data{6});
lon_sec    = double(data{7});

% Convert coordinates to decimal degrees
lat_decimal = -(lat_deg + lat_min / 60 + lat_sec / 3600);   % South is negative
lon_decimal = lon_deg + lon_min / 60 + lon_sec / 3600;      % East is positive

% Combine into XYZ matrix
geoid(geoid == -999) = NaN;
xyz = [lon_decimal, lat_decimal, geoid];

% Write output file
xyzFile ='AUSGeoid2020_20180201_win.xyz';
writematrix(xyz, xyzFile, 'Delimiter', ' ', 'FileType', 'text');

fprintf('✅ XYZ file written: %s\n', xyzFile);
