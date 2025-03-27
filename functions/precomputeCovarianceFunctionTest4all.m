function covariance_griddedInterpolant = precomputeCovarianceFunctionTest4all(covarianceType, BjerhammarRadius, ...
    widthDegree, resolutionDegree, coefficientA, coefficientB, maxOrder, minOrder)
    % This function computes the covariance between a signal at any point and outputs 
    % interpolatable functions for any distance and two radii, i.e., a full 3D spherical covariance function.
    %
    % Input: covarianceType =
    %              'cov_tt' = disturbing potential covariance model
    %              'cov_gg' = gravity anomaly covariance model
    %              'cov_gt' = gravity anomaly and disturbing potential cross covariance model
    %              'cov_dgdg' = gravity gradient and disturbing potential cross covariance model
    %              'cov_dgg' = gravity gradient covariance model
    %              'cov_dgt' = gravity gradient and gravity anomaly cross covariance model
    %         BjerhammarRadius = the radius of the Bjerhammar sphere
    %         widthDegree = size of the precomputed covariance function in degrees
    %         resolutionDegree = resolution of the covariance function in degrees
    %         coefficientA = covariance coefficient from computeCovarianceFunctionParameters
    %         coefficientB = covariance coefficient from computeCovarianceFunctionParameters
    %         maxOrder = max Legendre polynomial of covariance function
    %         minOrder = min Legendre polynomial of covariance function
    % Output: 
    %         covariance_griddedInterpolant = griddedInterpolant of the covariance model
    %           
    % Example: see computeLSC
    %
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11.

    resolutionDegree = deg2rad(resolutionDegree);
    widthDegree = deg2rad(widthDegree);
    
    % range of radius for topography of 500 m 
    psidat_int = 0:resolutionDegree:widthDegree;
    r1i = (BjerhammarRadius - 500:75:BjerhammarRadius + 500)';
    r2i = (BjerhammarRadius - 500:75:BjerhammarRadius + 500)';
    [r1im, r2im] = meshgrid(r1i, r2i);
    s = BjerhammarRadius^2 ./ (r1im .* r2im);
    % Feathestone numerical recepie, stable for high degree 
    pmm = ones(size(psidat_int));
    pmmp1 = cos(psidat_int) .* pmm;
    plgndr = pmmp1;
    
    covarianceModel = permute(repmat(plgndr, length(s(1, :)), 1, length(s(1, :))), [1, 3, 2]) * 0;

    CPsi = cos(psidat_int);

    for n = 2:maxOrder 
        
        pll = (CPsi .* (2 * n - 1) .* pmmp1 - (n - 1) .* pmm) ./ n;
        pmm = pmmp1;
        pmmp1 = pll;
        plgndr = pll;

        if n > minOrder   
            % Gravity anomaly degree variance model    
            degreeVariance = coefficientA * ((n - 1) / ((n - 2) * (n + coefficientB)));

            % Compute covariance based on the type
            switch covarianceType
                case 'cov_tt'
                    % Disturbing potential covariance model
                    covarianceModel = covarianceModel + degreeVariance * (BjerhammarRadius^2) / ((n - 1)^2) ...
                        * (repmat(s.^(n + 1), 1, 1, length(plgndr))) ...
                        .* permute(repmat(plgndr, length(s(1, :)), 1, length(s(1, :))), [1, 3, 2]);
                case 'cov_gg'
                    % Gravity anomaly covariance model
                    covarianceModel = covarianceModel + degreeVariance ...
                        * (repmat(s.^(n + 2), 1, 1, length(plgndr))) ...
                        .* permute(repmat(plgndr, length(s(1, :)), 1, length(s(1, :))), [1, 3, 2]);
                case 'cov_gt'
                    % Gravity anomaly and disturbing potential cross covariance model
                    covarianceModel = covarianceModel + degreeVariance / (n - 1) ...
                        * (repmat(s.^(n + 2), 1, 1, length(plgndr)) .* repmat(r2im, 1, 1, length(plgndr))) ...
                        .* permute(repmat(plgndr, length(s(1, :)), 1, length(s(1, :))), [1, 3, 2]);
                case 'cov_dgdg'
                    % Gravity gradient and disturbing potential cross covariance model
                    covarianceModel = covarianceModel + degreeVariance / (BjerhammarRadius^2) * ((n + 2)^2) ...
                        * (repmat(s.^(n + 3), 1, 1, length(plgndr))) ...
                        .* permute(repmat(plgndr, length(s(1, :)), 1, length(s(1, :))), [1, 3, 2]);
                case 'cov_dgg'
                    % Gravity gradient covariance model
                    covarianceModel = covarianceModel + degreeVariance / (n - 1) * ((n + 1) * (n + 2)) ...
                        * (repmat(s.^(n + 2), 1, 1, length(plgndr))) .* repmat(1 ./ (r1im), 1, 1, length(plgndr))...
                        .* permute(repmat(plgndr, length(s(1, :)), 1, length(s(1, :))), [1, 3, 2]);
                case 'cov_dgt'
                    % Gravity gradient and gravity anomaly cross covariance model
                    covarianceModel = covarianceModel + degreeVariance * (BjerhammarRadius^2) / ((n - 1)^2) * ((n + 1)^2) ...
                        * (repmat(s.^(n + 1), 1 , 1, length(plgndr)) .* repmat(1 ./ ((r1im.^2)), 1, 1, length(plgndr)))...
                        .* permute(repmat(plgndr, length(s(1, :)), 1, length(s(1, :))), [1, 3, 2]);
               otherwise
                   disp('Unexpected covariance type. No covariance model created.')
            end
        end
    end

    [r1imn, r2imn, psidat_intmn] = meshgrid(r1i, r2i, psidat_int);

    P = [2 1 3];
    r1imn = permute(r1imn, P);
    r2imn = permute(r2imn, P);
    psidat_intmn = permute(psidat_intmn, P);
    CCov_tt_intn = permute(covarianceModel, P);
   
    % make griddedInterpolant for each model
    covariance_griddedInterpolant = griddedInterpolant(double(r1imn), double(r2imn), double(psidat_intmn), CCov_tt_intn);

end
