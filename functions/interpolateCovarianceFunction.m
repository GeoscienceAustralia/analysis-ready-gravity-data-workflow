function covarianceMatrix = interpolateCovarianceFunction(longitude1, latitude1, height1, longitude2, latitude2, height2, covariance_griddedInterpolant)
    % interpolateCovarianceFunction computes the covariance between a signal at any point by
    % interpolating the 3D covariance function.
    %
    % Input:  longitude1, latitude1, height1 = vectors of coordinates of the first points
    %         longitude2, latitude2, height2 = vectors of coordinates of the second points
    %         covariance_griddedInterpolant = gridded interpolant for the covariance function
    %
    % Output: covarianceMatrix = matrix of interpolated covariances
    %           
    % Example: see computeLSC
    %
    % Main functions
    % - haversine
    %
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11.

    % Initialize covariance matrix
    covarianceMatrix = zeros(length(longitude1), length(longitude2));

    % Convert degrees to radians
    latitude1 = deg2rad (latitude1);
    latitude2 = deg2rad (latitude2);
    longitude1 = deg2rad (longitude1);
    longitude2 = deg2rad (longitude2);
    
    for k = 1:length(longitude1)
        haversineDistance = haversine(latitude1(k), longitude1(k), latitude2(:), longitude2(:));
        % Interpolate covariance function
        covarianceMatrix(k, :) = covariance_griddedInterpolant(double(height1(k)) * ones(size(height2')), double(height2'), double(haversineDistance'));
    end
     % Set NaN values to zero
    covarianceMatrix(isnan(covarianceMatrix)) = 0;
end
