function [coeff, z_detrended] = remove_tilted_plane_geo(lon_deg, lat_deg, z)
%REMOVE_TILTED_PLANE_GEO Remove best-fit tilted plane from z in geographic coords.
%
%   [coeff, z_detrended] = remove_tilted_plane_geo(lon_deg, lat_deg, z)
%
% Inputs:
%   lon_deg, lat_deg, z : arrays of same size (vectors or matrices)
%
% Outputs:
%   coeff       : [a b c] plane coefficients such that
%                 plane = a*(lon*cosd(lat)) + b*lat + c
%   z_detrended : z - plane (NaNs preserved where inputs are non-finite)

    % Ensure double precision (similar to np.asarray(..., float))
    lon = double(lon_deg);
    lat = double(lat_deg);
    z   = double(z);

    % Mask finite values only (equivalent to np.isfinite on all arrays)
    mask = isfinite(lon) & isfinite(lat) & isfinite(z);

    % Scale longitude by cos(latitude) (lat in degrees)
    lon_scaled = lon .* cosd(lat);

    % Build design matrix G = [lon_scaled, lat, 1]
    G = [lon_scaled(mask), lat(mask), ones(nnz(mask), 1)];

    % Solve least squares G*m ≈ z
    m = G \ z(mask);

    a = m(1); b = m(2); c = m(3);
    coeff = [a, b, c];

    % Compute plane everywhere (including where mask is false)
    plane = a .* lon_scaled + b .* lat + c;

    % Detrend; keep NaNs where original z is non-finite
    z_detrended = z - plane;
    z_detrended(~isfinite(z)) = NaN;
end