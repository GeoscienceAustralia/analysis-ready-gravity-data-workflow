function C = mixed_gaussian_covariance(psi, c, L)
%MIXED_GAUSSIAN_COVARIANCE Mixed Gaussian covariance function.
%
%   C = mixed_gaussian_covariance(psi, c, L)
%
% Computes:
%   C(psi) = sum_i c(i) * exp( -0.5 * (psi / L(i))^2 )
%
% Inputs
%   psi : array (any shape), angular separation in radians
%   c   : vector (Kx1 or 1xK), positive covariance amplitudes
%   L   : vector (Kx1 or 1xK), positive correlation lengths (radians)
%
% Output
%   C   : array, same size as psi, covariance evaluated at psi
%
% Notes
%   - Fully vectorized over psi and mixture components i.
%   - Preserves the shape of psi in the output.

    % --- Type/shape normalization ---
    psi = double(psi);            % allow any shape
    c   = double(c(:));           % force column vector
    L   = double(L(:));           % force column vector

    % --- Validation ---
    if numel(c) ~= numel(L)
        error('mixed_gaussian_covariance:SizeMismatch', ...
              'c and L must have the same length.');
    end
    if any(~isfinite(c)) || any(~isfinite(L)) || any(~isfinite(psi(:)))
        error('mixed_gaussian_covariance:NonFinite', ...
              'psi, c, and L must be finite.');
    end
    if any(c <= 0)
        error('mixed_gaussian_covariance:NonPositiveC', ...
              'All c values must be > 0.');
    end
    if any(L <= 0)
        error('mixed_gaussian_covariance:NonPositiveL', ...
              'All L values must be > 0.');
    end

    % --- Vectorized computation ---
    % Make psi a column for broadcasting, then reshape back.
    psi_col = psi(:);                         % (N x 1)
    % Compute (psi/L)^2 for all N,K using implicit expansion:
    % psi_col is (N x 1), L' is (1 x K) => (N x K)
    X = (psi_col ./ (L.')).^2;                % (N x K)
    E = exp(-0.5 .* X);                       % (N x K)
    C_col = E * c;                            % (N x 1) weighted sum over K

    % Reshape to original psi shape
    C = reshape(C_col, size(psi));
end