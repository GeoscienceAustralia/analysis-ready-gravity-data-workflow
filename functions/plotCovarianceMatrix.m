function plotCovarianceMatrix(distanceRad, ...
    covariance, covarianceUnit, ...
    covarianceTitle, outputParameters)
% plotSphericalCovarianceFunction
% Plots empirical and fitted spherical covariance functions
%
% Inputs:
%   distanceRad       - Spherical distance in radians
%   covariance        - Empirical covariance values
%   fittedCovariance  - Fitted covariance values
%   covarianceUnit    - Y-axis label (LaTeX-compatible)
%   covarianceTitle   - Figure title (string)
%   outputParameters  - Output path or prefix for saved figure

    % Convert distance to degrees for plotting
    distanceDeg = rad2deg(distanceRad);

    % Create figure
    figure('Name','Covariance Function','NumberTitle','off');
    clf
    hold on
    box on
    grid on

    % Plot data
    plot(distanceDeg, covariance(end,:), 'k*')
    
    % Axes formatting
    xlabel('Spherical distance (degrees)')
    xlim([0 1.2])
    ylabel(covarianceUnit, 'Interpreter','latex')
    title(covarianceTitle, 'Interpreter','none')
    ax = gca;
    ax.YAxis.Exponent = 0;              % remove ×10^n
    ax.YAxis.TickLabelFormat = '%.3f';  % fixed decimal format
    
    drawnow

    % Make filename safe
    safeTitle = regexprep(covarianceTitle,'[^\w]','_');
    filename = fullfile(outputParameters, ...
        ['Covariance_', safeTitle, '.png']);
    % Save figure
    saveas(gcf, filename)


     figure
     imagesc(covariance)
     set(gca,'YDir','normal')
     colorbar
     colormap(jet)
     title(colorbar,covarianceUnit,'FontSize',10);
     title(covarianceTitle, 'Interpreter','none')
     % Save figure
     filename = fullfile(outputParameters, ...
        ['CovarianceMatrix', safeTitle, '.png']);
     saveas(gcf,filename)

end
