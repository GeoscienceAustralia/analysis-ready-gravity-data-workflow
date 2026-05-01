function plotSphericalCovarianceFunction(distanceRad, ...
    covariance, fittedCovariance, covarianceUnit, ...
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
    plot(distanceDeg, covariance, 'k*', 'DisplayName','Empirical data')
    plot(distanceDeg, fittedCovariance, 'r-', ...
        'LineWidth', 1.5, 'DisplayName','Fitted function')

    % Axes formatting
    xlabel('Spherical distance in degrees')
    xlim([0 max(distanceDeg)])
    ylabel('Covariance (m^2)')
    ylabel(covarianceUnit, 'Interpreter','latex')
    title(covarianceTitle, 'Interpreter','none')
    ax = gca;
    ax.YAxis.Exponent = 0;              % remove ×10^n
    ax.YAxis.TickLabelFormat = '%.3f';  % fixed decimal format
    legend('Location','best')

    drawnow

    % Make filename safe
    safeTitle = regexprep(covarianceTitle,'[^\w]','_');
    filename = fullfile(outputParameters, ...
        ['Covariance_', safeTitle, '.png']);

    % Save figure
    saveas(gcf, filename)

end
