function fig = plotGridedGeometriCorrection(Zq,LongDEM, LatDEM,climits, ...
        GRID_PARA, OUTPUT_PARA, Coastline,title)
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
    % Axis limits
    axisLimits.latMeanCosine = abs(cosd(mean([GRID_PARA.MINLAT, GRID_PARA.MAXLAT])));
    axisLimits.lonMinLimit   = GRID_PARA.MINLONG - GRID_PARA.buffer;
    axisLimits.lonMaxLimit   = GRID_PARA.MAXLONG + GRID_PARA.buffer;
    axisLimits.latMinLimit   = GRID_PARA.MINLAT  - GRID_PARA.buffer;
    axisLimits.latMaxLimit   = GRID_PARA.MAXLAT  + GRID_PARA.buffer;

    % --------------------------------------------------------------
    % Robust colour limits (±2σ)
    
    if isempty(climits)
        mu  = mean(Zq(:), 'omitnan');
        sig = std(Zq(:),  'omitnan');
    
        cmin = mu - 2*sig;
        cmax = mu + 2*sig;
    else
        cmin = climits(1);
        cmax = climits(2);
    end

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

    customizeMap(title,'m',Coastline,axisLimits)

    % --------------------------------------------------------------
    % Save output
    if ~exist(OUTPUT_PARA.plotsFolder, 'dir')
        mkdir(OUTPUT_PARA.plotsFolder)
    end

    outName = fullfile( ...
        OUTPUT_PARA.plotsFolder, ...
        [title,'.png']);

    exportgraphics(fig, outName, 'Resolution', 300)

end