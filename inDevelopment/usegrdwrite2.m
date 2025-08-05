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


this is how i wrote the data in matlab:


% Combine columns into a single matrix
    data = [resamplegeoid_col, Lat_out_col_DM, zero_columns, Long_out_col_DM, zero_columns, zero_columns, zero_columns];
    
    % Open the file for writing
    fileID = fopen([OUTPUT_PARA.Grids_name, 'AGQG_', formattedDate, '.dat'], 'w');
    
    % Write header
    fprintf(fileID, '../data/AGQG_%s                           www.ga.gov.au\n', formattedDate);
    
    % Write data in the specified format
    %for i = 1:8872
    for i = 1:size(data,1)
        fprintf(fileID, 'GEO  %7.3f S%2d %2.f  %.3f E%3d %2.f  %.3f      %2.2f      %2.2f\n', ...
                data(i, 1), data(i, 2), data(i, 3), data(i, 4), data(i, 5), data(i, 6), data(i, 7), data(i, 8), data(i, 9));
    end
    

i want to read this file and write into a xyz file with three columns longtud in decimal, latitude in decimal and resamplegeoid_col.
