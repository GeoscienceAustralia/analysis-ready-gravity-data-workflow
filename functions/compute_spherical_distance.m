function psi = compute_spherical_distance(lat1, lon1, lat2, lon2)
%COMPUTE_SPHERICAL_DISTANCE
% Compute angular (great-circle) distance between two point sets on a sphere.
%
% Inputs
%   lat1, lon1 : vectors (n1 x 1 or 1 x n1), degrees
%   lat2, lon2 : vectors (n2 x 1 or 1 x n2), degrees
%
% Output
%   psi : (n1 x n2) matrix of angular separations [radians]
%
% Notes
%   - Uses the haversine formula
%   - Fully vectorised (no loops)
%   - Output is in radians

    % Ensure column vectors
    lat1 = lat1(:);
    lon1 = lon1(:);
    lat2 = lat2(:);
    lon2 = lon2(:);

    % Convert degrees to radians
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    % Expand dimensions for broadcasting:
    % lat1, lon1 -> (n1 x 1)
    % lat2, lon2 -> (1 x n2)
    dlat = lat2.' - lat1;   % (n1 x n2)
    dlon = lon2.' - lon1;   % (n1 x n2)

    % Haversine formula
    a = sin(dlat ./ 2).^2 + ...
        cos(lat1) * cos(lat2.') .* sin(dlon ./ 2).^2;

    % Numerical safety (optional but good practice)
    a = min(max(a, 0), 1);

    psi = 2 .* asin(sqrt(a));
end