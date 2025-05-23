function covarianceParameters = computeCovarianceFunctionParameters(outputParameters, covarianceParameters, gravityData, gravityType, block_counter)
    % Compute spherical empirical covariance and fit covariance function to empirical covariance.
    % No gradient data are used here at all.
    %
    % Input:  gridParameters = grid/tiling parameters such as extent (MINLONG, MAXLONG, MINLAT, MAXLAT), buffer, STEP 
    %         covarianceParameters = covariance function parameters
    %         gravityData = [longitude latitude orthometric_height gravity_anomaly uncertainty(std) flag(1:terrestrial 2:airborne 3:gradiometry)]
    %         gravityType = character, e.g. Faye, RTM
    %         block_counter = block number 
    % Output: covarianceParameters = updated covariance parameters
    %           
    % Example: see computeLSC
    %
    % Main functions
    % - computeSphericalEmpiricalCovariance 
    % - fitEmpiricalCovariance
    % - plotCovarianceFunction
    %
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11.

    if covarianceParameters.Airbornedataonly
        gravityData = gravityData(gravityData(:,6) == 2, :);
    end

    residuals = gravityData(:,4); 
    residuals = residuals - mean(residuals);

    % Calculate empirical covariances
    covarianceInfo = computeSphericalEmpiricalCovariance(gravityData(:,1), gravityData(:,2), residuals, covarianceParameters.Compute_Empircal_COV_Dec);

    % Fit empirical covariances
    cont = 0;
    if strcmp(covarianceParameters.Fit_Empircal_COV, 'man')
        while cont == 0
            [covarianceParameters.A, covarianceParameters.B, fittedCovariance] = fitEmpiricalCovariance(covarianceInfo(:,1), covarianceInfo(:,2), covarianceParameters.N, covarianceParameters.M);
            
            figure(17)
            clf
            hold on
            plot(rad2deg ( covarianceInfo(:,1) ), covarianceInfo(:,2), '*')
            plot(rad2deg ( covarianceInfo(:,1) ), covarianceParameters.A * fittedCovariance, '-')
            drawnow
            legend('Empirical data', 'Fitted function')
            xlabel('Spherical distance in degrees')
            title('Gravity auto-covariance')
            
            disp('Done')
            COVInput = input('Try different [N, M]? 1/0 = y/n: ');
            if COVInput == 1
                covarianceParameters.N = input('N? ');
                covarianceParameters.M = input('M? ');
            else
                cont = 1;
            end
        end
        saveas(gcf, 'covariance.png')    
    elseif strcmp(covarianceParameters.Fit_Empircal_COV, 'auto')
        disp(['Fitting covariance function for ',gravityType])
        bestEstimate = [1e6, 1e6];
        
        for N = covarianceParameters.FitEmpiricalCOVNSearch(1):covarianceParameters.FitEmpiricalCOVNSearch(2):covarianceParameters.FitEmpiricalCOVNSearch(3)
            for M = covarianceParameters.FitEmpiricalCOVMSearch(1):covarianceParameters.FitEmpiricalCOVMSearch(2):covarianceParameters.FitEmpiricalCOVMSearch(3)
                [A, B, fittedCovariance] = fitEmpiricalCovariance(covarianceInfo(:,1), covarianceInfo(:,2), N, M);

                residualCov = covarianceInfo(:,2) - interp1(covarianceInfo(:,1), A * fittedCovariance, covarianceInfo(:,1));
                EST = residualCov;
                
                if sum(abs(EST)) < sum(abs(bestEstimate))
                    bestEstimate = EST;
                    covarianceParameters.N = N;
                    covarianceParameters.M = M;
                    covarianceParameters.A = A;
                    covarianceParameters.B = B;
                end
            end
        end
    else
        [covarianceParameters.A, covarianceParameters.B, fittedCovariance] = fitEmpiricalCovariance(covarianceInfo(:,1), covarianceInfo(:,2), covarianceParameters.N, covarianceParameters.M);
    end

    if covarianceParameters.COVPlot
    %plotCustomScatter(gravityData(:,1),gravityData(:,2),residuals,gridParameters,1,gravityType,'mGal','Outputs/plots/')
    plotCovarianceFunction(gravityData, residuals, outputParameters, gravityType, covarianceInfo, covarianceParameters, fittedCovariance,block_counter);
    end
end
