function plotSphericalCovarianceFunction(distanceDegree, covariance, fittedCovariance,covarianceUnit,covarianceTitle, outputParameters)

    % Plot covariance function
    
    figure('Name','computeCovarianceFunction','NumberTitle','off');
    clf
    hold on
    plot(rad2deg ( distanceDegree ), covariance, '*')
    plot(rad2deg ( distanceDegree ), fittedCovariance, '-')
    xlim([0 1]);
    drawnow
    legend('Empirical data', 'Fitted function')
    xlabel('Spherical distance in degrees')
    ylabel(covarianceUnit, 'interpreter', 'latex')
    title(covarianceTitle)   
    saveas(gcf,[outputParameters,covarianceTitle,'.png'])

end

