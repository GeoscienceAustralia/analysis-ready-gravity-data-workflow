function covarianceInfo = computeSphericalEmpiricalCovariance(Long, Lat, Gravity, DecimationFactor)
    % Computes the empirical covariance for spherical distance
    %
    % Input:  Long = longitude vector
    %         Lat = latitude vector
    %         Gravity = gravity anomaly vector
    %         DecimationFactor = Decimation factor for empirical covariance estimates. e.g. 1 is no decimation, 2 drops 50% of the data etc.
    %
    % Output: covarianceInfo= [covdist',covariance];
    %
    % Example: see computeCovarianceFunctionParameters
    %
    % Main functions
    % - custom_grpstats
    % - haversine
    % Other functions
    % -
    % -
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11.

    if length(Gravity) < 120
    DecimationFactor = 2;
    end
    
    maxSphericalDistance = deg2rad(1);
    ds = deg2rad(2/60);
    covariance = zeros(round(maxSphericalDistance/ds) + 2, 1);
    ncov = zeros(round(maxSphericalDistance/ds) + 2, 1);
    covdist = (0:length(covariance)-1) * ds;
    
    % Some decimation
    % Fix the seed for reproducibility
    rng(42);  % You can use any integer as the seed
    randNumbers=randn(size(Lat));
    [~, ind] = sort(randNumbers);
    ind(1:end - round(length(Lat)/DecimationFactor)) = [];
    Dataset.Lat = Lat(ind);
    Dataset.Long = Long(ind);
    gravityAnomalies = Gravity(ind);

    % Convert degrees to radians
    Dataset.Lat = deg2rad (Dataset.Lat);
    Dataset.Long = deg2rad (Dataset.Long);
    
    % Computes the empirical covariance
    for i = 1:length(gravityAnomalies)
        haversineDist = haversine(Dataset.Lat(i), Dataset.Long(i), Dataset.Lat(:), Dataset.Long(:));
        ir = round((haversineDist/ds)) + 1;
        Tooadd = gravityAnomalies(i) * gravityAnomalies(:);
        Tooadd(haversineDist > maxSphericalDistance) = 0;
        [sumir, count] = custom_grpstats(Tooadd, ir, {@sum, @numel});
        covariance(unique(ir(ir <= length(covariance)))) = covariance(unique(ir(ir <= length(covariance)))) + sumir(unique(ir(ir <= length(covariance))));
        ncov(unique(ir(ir <= length(ncov)))) = ncov(unique(ir(ir <= length(ncov)))) + count(unique(ir(ir <= length(ncov))));
    end
    
    for i = 1:length(covariance)
        if ncov(i) ~= 0
            covariance(i) = covariance(i) / ncov(i);
        end
    end
    
    covarianceInfo = [covdist', covariance];
    
end
