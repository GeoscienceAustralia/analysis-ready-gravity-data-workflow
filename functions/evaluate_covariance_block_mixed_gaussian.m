function C_block = evaluate_covariance_block_mixed_gaussian( ...
    lat1, lon1, lat2, lon2, interp_C)
%EVALUATE_COVARIANCE_BLOCK_MIXED_GAUSSIAN
% Evaluate covariance block between two point sets using
% a psi-based covariance interpolator.
%
% Inputs
%   lat1, lon1 : vectors (n1 x 1) or (1 x n1)
%   lat2, lon2 : vectors (n2 x 1) or (1 x n2)
%   interp_C   : griddedInterpolant or function handle
%
% Output
%   C_block    : (n1 x n2) covariance matrix

    % Ensure column vectors
    lat1 = lat1(:);
    lon1 = lon1(:);
    lat2 = lat2(:);
    lon2 = lon2(:);

    % Compute spherical distance matrix (n1 x n2)
    psi = compute_spherical_distance(lat1, lon1, lat2, lon2);

    % Flatten, evaluate interpolator, and reshape
    psi_flat = psi(:);
    C_flat   = interp_C(psi_flat);
    C_block  = reshape(C_flat, numel(lat1), numel(lat2));
end