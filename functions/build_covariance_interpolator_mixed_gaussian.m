function interp = build_covariance_interpolator_mixed_gaussian(psi_vals, c, L)
%BUILD_COVARIANCE_INTERPOLATOR_MIXED_GAUSSIAN
% Build a 1D interpolator in psi for the mixed Gaussian covariance.
%
% Equivalent intent to Python/SciPy:
%   psi_vals = np.asarray(psi_vals, float)
%   C_grid   = build_covariance_grid_mixed_gaussian(psi_vals, c, L)
%   interp   = RegularGridInterpolator((psi_vals,), C_grid,
%                                     bounds_error=False, fill_value=None)
%
% Inputs:
%   psi_vals : vector of grid points (radians)
%   c, L     : mixture parameters (vectors, same length)
%
% Output:
%   interp   : griddedInterpolant object; evaluate as interp(psi_query)
%              Extrapolates linearly outside the grid.

    % Normalize inputs
    psi = double(psi_vals(:));  % column vector
    if numel(psi) < 2
        error('psi_vals must contain at least 2 points for interpolation.');
    end

    % Sort grid (SciPy expects monotonic grid; MATLAB also benefits)
    [psi_sorted, idx] = sort(psi);

    % Build covariance grid on sorted psi
    C_grid = build_covariance_grid_mixed_gaussian(psi_sorted, c, L);
    C_grid = double(C_grid(:));

    if numel(C_grid) ~= numel(psi_sorted)
        error('Covariance grid must have the same number of points as psi_vals.');
    end

    % Create interpolant:
    % - 'linear' interpolation (typical default)
    % - 'linear' extrapolation to mirror fill_value=None (extrapolate)
    interp = griddedInterpolant(psi_sorted, C_grid, 'linear', 'linear');
end