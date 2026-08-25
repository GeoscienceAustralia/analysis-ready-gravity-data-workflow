function C_grid = build_covariance_grid_mixed_gaussian(psi_vals, c, L)
%BUILD_COVARIANCE_GRID_MIXED_GAUSSIAN
% Build a 1D covariance grid over psi only.
%
% Returns:
%   C_grid : vector, same length as psi_vals

    psi = double(psi_vals(:));
    C_grid = mixed_gaussian_covariance(psi, c, L);  % you already translated this
end
