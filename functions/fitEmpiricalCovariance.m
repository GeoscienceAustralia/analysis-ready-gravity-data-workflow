function [bestFitCoeff_A, bestFitCoeff_B, fittedCovariance] = fitEmpiricalCovariance(sphericalDistance, empiricalCovariance, maxOrder, minOrder)
    % Fits empirical covariance to Legendre polynomial expansion.
    %
    % Reference: Equation(68) in Tscherning 1974, page 30
    % 
    % Input:  sphericalDistance = vector of spherical distances
    %         empiricalCovariance = vector of empirical covariance
    %         maxOrder = maximum order for Legendre polynomials
    %         minOrder = minimum order for Legendre polynomials
    % 
    % Output: bestFitCoeff_A = best-fit coefficient A
    %         bestFitCoeff_B = best-fit coefficient B
    %         fittedCovariance = vector of fitted covariance
    %           
    % Example: see computeCovarianceFunction
    %
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11.

    constants                                       % load constants
    % Constants
    %EarthRadius = 6371000; % is this from Tscherning 1974 abstract?
    maxIterations = 250;

    % Set coefficients
    bestFitCoeff_A = 1;
    bestFitCoeff_B = 0;

    for iteration = 1:maxIterations
        % Compute the function
        scale_factor = 1;

        pmm = ones(size(sphericalDistance));
        pmmp1 = cos(sphericalDistance) .* pmm;
        plgndr = pmmp1;

        C2_da = 0 * (scale_factor .^ 2) .* plgndr;
        C2_db = 0 * (scale_factor .^ 2) .* plgndr;
        computedCovariance = 0 * (scale_factor .^ 2) .* plgndr;
        fittedCovariance = 0 * (scale_factor .^ 2) .* plgndr;

        for index = 2:maxOrder
            
            pll = (cos(sphericalDistance) .* (2 * index - 1) .* pmmp1 - (index - 1) .* pmm) / index;
            pmm = pmmp1;
            pmmp1 = pll;
            plgndr = pll;

            if index >= minOrder
                fittedCovariance = fittedCovariance + bestFitCoeff_A * ((index - 1) / ((index - 2) * (index + bestFitCoeff_B))) * (scale_factor .^ (index + 2)) .* plgndr;    
                C2_da = C2_da + ((index - 1) / ((index - 2) * (index + bestFitCoeff_B))) * (scale_factor .^ (index + 2)) .* plgndr;
                C2_db = C2_db - bestFitCoeff_A * ((index - 1) / ((index - 2) * ((index + bestFitCoeff_B) .^ 2))) * (scale_factor .^ (index + 2)) .* plgndr;
                computedCovariance = computedCovariance + (EarthRadius ^ 2) / ((index - 1) ^ 2) * ((index - 1) / ((index - 2) * (index + bestFitCoeff_B))) * (scale_factor .^ (index + 2)) .* plgndr;
            end
        end

        % Compute best fitted coefficients
        if iteration == 1
            bestFitCoeff_A = C2_da \ empiricalCovariance;
        else
            Sol = [C2_da, C2_db] \ (empiricalCovariance - fittedCovariance);
            bestFitCoeff_A = abs(bestFitCoeff_A + 0.1 * Sol(1));
            bestFitCoeff_B = bestFitCoeff_B + 0.2 * Sol(2);
        end
    end
    fittedCovariance = fittedCovariance / bestFitCoeff_A;
end