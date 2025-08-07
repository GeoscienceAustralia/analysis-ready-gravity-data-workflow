function xyz = readDatFile(inputFile)
    % Check file exists
    fid = fopen(inputFile, 'r');
    if fid == -1
        error('Failed to open file: %s', inputFile);
    end

    % Read data (adjust format if necessary)
    data = textscan(fid, 'GEO %f S%f %f %f E%f %f %f %f %f', 'HeaderLines', 1);
    fclose(fid);

    % Extract values
    geoid     = double(data{1});
    lat_deg   = double(data{2});
    lat_min   = double(data{3});
    lat_sec   = double(data{4});
    lon_deg   = double(data{5});
    lon_min   = double(data{6});
    lon_sec   = double(data{7});

    % Convert coordinates to decimal degrees
    lat_decimal = -(lat_deg + lat_min / 60 + lat_sec / 3600);   % South
    lon_decimal =  (lon_deg + lon_min / 60 + lon_sec / 3600);   % East

    % Handle no-data
    geoid(geoid == -999) = NaN;

    % Combine into output
    xyz = [lon_decimal, lat_decimal, geoid];

    % Plot
    figure;
    scatter(lon_decimal, lat_decimal, 1, geoid, 'filled');
    colorbar;
    colormap(parula);
    xlabel('Longitude');
    ylabel('Latitude');
    title('Geoid in meters');
end

