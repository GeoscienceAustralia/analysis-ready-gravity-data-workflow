function fig = plotInterpolatedGeometricCorrection( ...
        Lev, geomGravGeoidDiffDetrended, ...
        LongDEM, LatDEM, ...
        GRID_PARA, OUTPUT_PARA, Coastline)
% plotInterpolatedGeometricCorrection
% Interpolates detrended gravity/geoid residuals and plots a mosaic map
%
% Inputs:
%   Lev                         [Nx2]  lon, lat
%   geomGravGeoidDiffDetrended  [Nx1]  residual values (m)
%   LongDEM, LatDEM             DEM grids (meshgrid-style)
%   GRID_PARA                   struct with MIN/MAX LAT/LONG and buffer
%   OUTPUT_PARA                 struct with plotsFolder
%   Coastline                   coastline struct for customizeMap
%
% Output:
%   fig                         figure handle

    % --------------------------------------------------------------
    % Extract coordinates
    lon = Lev(:,1);
    lat = Lev(:,2);
    z   = geomGravGeoidDiffDetrended;

    % --------------------------------------------------------------
    % Remove invalid values
    idx = isfinite(lon) & isfinite(lat) & isfinite(z);
    lon = lon(idx);
    lat = lat(idx);
    z   = z(idx);

    % --------------------------------------------------------------
    % Create interpolant
    F = scatteredInterpolant( ...
        lon, lat, z, ...
        'natural', ...   % interpolation
        'none');         % no extrapolation

    % Evaluate on DEM grid
    Zq = F(LongDEM, LatDEM);

    % --------------------------------------------------------------
    % Axis limits
    axisLimits.latMeanCosine = abs(cosd(mean([GRID_PARA.MINLAT, GRID_PARA.MAXLAT])));
    axisLimits.lonMinLimit   = GRID_PARA.MINLONG - GRID_PARA.buffer;
    axisLimits.lonMaxLimit   = GRID_PARA.MAXLONG + GRID_PARA.buffer;
    axisLimits.latMinLimit   = GRID_PARA.MINLAT  - GRID_PARA.buffer;
    axisLimits.latMaxLimit   = GRID_PARA.MAXLAT  + GRID_PARA.buffer;

    % --------------------------------------------------------------
    % Robust colour limits (±2σ)
    mu  = mean(Zq(:), 'omitnan');
    sig = std(Zq(:),  'omitnan');

    cmin = mu - 2*sig;
    cmax = mu + 2*sig;

    % --------------------------------------------------------------
    % Plot
    fig = figure( ...
        'Name','Geometric correction', ...
        'NumberTitle','off', ...
        'Color','w');

    clf
    hold on

    imagesc(LongDEM(1,:), LatDEM(:,1), Zq)
    set(gca,'YDir','normal')
    colormap(jet)
    caxis([cmin cmax])

    customizeMap('Geometric correction','m',Coastline,axisLimits)

    % --------------------------------------------------------------
    % Save output
    if ~exist(OUTPUT_PARA.plotsFolder, 'dir')
        mkdir(OUTPUT_PARA.plotsFolder)
    end

    outName = fullfile( ...
        OUTPUT_PARA.plotsFolder, ...
        'MosaicTiles_Geometric_correction.png');

    exportgraphics(fig, outName, 'Resolution', 300)

end