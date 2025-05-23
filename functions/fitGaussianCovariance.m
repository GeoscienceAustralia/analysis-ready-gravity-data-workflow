function [sigma2,bestFitCoeff,fittedCovariance] = fitGaussianCovariance(sphericalDistance,empiricalCovariance)
    % Fits empirical covariance to Gaussian function.
    %
    % Reference: Brown et. al (2018) page 1458 Fig.2
    % 
    % Input:  sphericalDistance = vector of spherical distances
    %         empiricalCovariance = vector of empirical covariance
    %
    % Output: bestFitCoeff_A = best-fit coefficient A
    %         bestFitCoeff_B = best-fit coefficient B
    %         fittedCovariance = vector of fitted covariance
    %           
    % Example: see computeCovarianceFunction
    %
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2024-07.
    
    % Set coefficients
    sigma2 = 1;
    bestFitCoeff = 0.001;

    for k=1:250
        % Compute the function and the derivatives
        fittedCovariance = sigma2*exp(-(sphericalDistance.^2)/(2*bestFitCoeff^2));    
        C2_dsigma2 = exp(-(sphericalDistance.^2)/(2*bestFitCoeff^2));
        C2_dl = (sphericalDistance.^2)./(bestFitCoeff^3).*sigma2.*exp(-(sphericalDistance.^2)/(2*bestFitCoeff^2));
        % Compute best fitted coefficients
        if k == 1
            sigma2 = C2_dsigma2\empiricalCovariance;
            else
            solutionVector = [C2_dsigma2(1:5),C2_dl(1:5)]\(empiricalCovariance(1:5)-fittedCovariance(1:5));
            sigma2 = abs(sigma2+0.1*solutionVector(1));
            bestFitCoeff = bestFitCoeff+0.2*solutionVector(2);
        end
    end
    
end