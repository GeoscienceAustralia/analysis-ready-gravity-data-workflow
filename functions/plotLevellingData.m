function fig = plotLevellingData(Lon, Lat, Lev, quantityName, plotsFolder)
% plotLevellingData  Plot GPS levelling data on a lon/lat scatter map
%
% Inputs:
%   Lon           [Nx1] longitude (deg)
%   Lat           [Nx1] latitude (deg)
%   Lev           [Nx1] levelling values (m)
%   quantityName  string/char, plot title
%   plotsFolder   string/char (optional), output folder
%
% Output:
%   fig           figure handle

    arguments
        Lon (:,1) double
        Lat (:,1) double
        Lev (:,1) double
        quantityName {mustBeText}
        plotsFolder {mustBeText} = ""
    end

    % --- Robust colour limits (match your Python logic) ---
    mu  = mean(Lev, 'omitnan');
    sig = std(Lev,  'omitnan');

    cmin = mu - 2*sig;
    cmax = mu + 2*sig;

    % --- Create figure ---
    fig = figure( ...
        'Name', 'GPS Levelling', ...
        'NumberTitle', 'off', ...
        'Color', 'w');

    clf(fig)
    hold on

    % --- Scatter plot ---
    scatter(Lon, Lat, 10, Lev, 'filled')

    % --- Axes & appearance ---
    axis equal tight
    grid on
    colormap(jet)
    caxis([cmin cmax])
    set(gca, 'FontSize', 11)

    % --- Colorbar ---
    cb = colorbar;
    cb.Label.String  = 'm';
    cb.Label.FontSize = 11;

    % --- Labels & title ---
    xlabel('Longitude')
    ylabel('Latitude')
    title(quantityName, 'Interpreter','none')

    % --- Save outputs if requested ---
    if strlength(plotsFolder) > 0
        if ~exist(plotsFolder, 'dir')
            mkdir(plotsFolder)
        end

        % Make filename safe
        safeName = regexprep(quantityName, '[^a-zA-Z0-9_-]', '_');
        baseName = fullfile(plotsFolder, "levelling_" + safeName);

        savefig(fig, baseName + ".fig")
        exportgraphics(fig, baseName + ".png", 'Resolution', 300)
    end

end