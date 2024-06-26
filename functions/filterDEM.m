function [filtered_dem] = filterDEM(dem, harmonic_degree, lat_extent, lon_extent)
    % filterDEM provides a filtered DEM by Fourier transform and Low-pass filter mask

    % Input:  dem = 2D array representing the Digital Elevation Model
    %         harmonic_degree = the degree of the spherical harmonic model used as a reference
    %         lat_extent = [min_lat, max_lat] - Latitude extent of the DEM
    %         lon_extent = [min_lon, max_lon] - Longitude extent of the DEM
    %
    % 
    % Output:  filtered_dem = 2D array representing the filtered DEM
       
    % Example: see RTM_CORRECTIONS

    % Main functions
    % - 
    % Other functions
    % -   
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11.

    % Fourier transform of the DEM
    dem_fft = fft2(dem);

    % Create a frequency domain mask based on the harmonic degree
    [rows, cols] = size(dem);
    
    % Create frequency grids that are adaptable for even and odd sizes
    column_shifts = ifftshift(-ceil((cols-1)/2):floor(cols/2));
    row_shifts = ifftshift(-ceil((rows-1)/2):floor(rows/2));
    [freq_U, freq_V] = meshgrid(column_shifts, row_shifts);
    
    radius = sqrt(freq_U.^2 + freq_V.^2);  % Frequency magnitude
    
    % Compute the actual DEM dimension in harmonic degrees
    dem_dim_deg_lat = abs(lat_extent(2) - lat_extent(1));
    dem_dim_deg_lon = abs(lon_extent(2) - lon_extent(1));
    
    max_dem_dim = max([dem_dim_deg_lat, dem_dim_deg_lon]);
    
    % Adjust the cutoff based on DEM dimensions and desired harmonic degree
    scaled_radius = radius * (180 / max_dem_dim);
    cutoff = (harmonic_degree / 180) * 180;

    % Low-pass filter mask
    mask = scaled_radius <= cutoff;

    % Apply the mask to the Fourier-transformed DEM
    dem_fft_filtered = dem_fft .* mask;

    % Inverse Fourier transform to get the filtered DEM
    filtered_dem = real(ifft2(dem_fft_filtered));

end
