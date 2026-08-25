function [z_restored, plane_grid] = restore_tilted_plane_geo( ...
    lon_grid_deg, lat_grid_deg, z_detrended_grid, plane_params)
%RESTORE_TILTED_PLANE_GEO Restore a previously removed tilted plane on a grid.
%
% Inputs:
%   lon_grid_deg        : grid of longitudes (degrees)
%   lat_grid_deg        : grid of latitudes  (degrees)
%   z_detrended_grid   : detrended grid values
%   plane_params       : [a b c] from remove_tilted_plane_geo
%
% Outputs:
%   z_restored : grid with tilted plane restored
%   plane_grid : evaluated tilted plane on the grid

    % Ensure double precision (like np.asarray(..., float))
    lon_grid = double(lon_grid_deg);
    lat_grid = double(lat_grid_deg);
    z_detrended_grid = double(z_detrended_grid);

    % Unpack plane parameters
    a = plane_params(1);
    b = plane_params(2);
    c = plane_params(3);

    % Scale longitude by cos(latitude)
    lon_scaled_grid = lon_grid .* cosd(lat_grid);

    % Evaluate plane on the grid
    plane_grid = a .* lon_scaled_grid + b .* lat_grid + c;

    % Restore the plane
    z_restored = z_detrended_grid + plane_grid;
end