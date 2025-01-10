function DisplayAreaStatistics(Coastline,GRID_PARA,LongDEM,LatDEM,Grid_res_geoid_w, ...
    Grid_res_geoid_err_w,OUTPUT_PARA)

    % common variables for plotting
    axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
    axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
    axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
    axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
    axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;

    if ~isempty(OUTPUT_PARA.polygonLon) && ~isempty(OUTPUT_PARA.polygonLat)

       areaIn=inpolygon(LongDEM,LatDEM,OUTPUT_PARA.polygonLon,OUTPUT_PARA.polygonLat);
       Grid_res_geoid_err_w(areaIn==0)=nan;
       Grid_res_geoid_w(areaIn==0)=nan;
    
        figure('Name','MosaicTiles','NumberTitle','off'); 
        clf
        hold on
        imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid_err_w)
        hold on;
        plot(OUTPUT_PARA.polygonLon, OUTPUT_PARA.polygonLat, 'magenta-', 'LineWidth', 2); % Plot the polygon
        customizeMap('residualGeoidSigmaError','m',Coastline,axisLimits)
        %caxis([0 0.05])
        caxis([0 0.02])
        saveas(gcf,[OUTPUT_PARA.plotsFolder,'MosaicTiles','residualGeoidSigmaErrorFocused','.png'])
    
        figure('Name','MosaicTiles','NumberTitle','off'); 
        clf
        hold on
        imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid_w)
        hold on;
        plot(OUTPUT_PARA.polygonLon, OUTPUT_PARA.polygonLat, 'magenta-', 'LineWidth', 2); % Plot the polygon
        customizeMap('residualGeoid','m',Coastline,axisLimits)
        caxis([-0.5 0.5])
        saveas(gcf,[OUTPUT_PARA.plotsFolder,'MosaicTiles','residualGeoidFocused','.png'])

    end

    disp('Display statistics')

    Grid_res_geoid_err_w(isnan(Grid_res_geoid_err_w)) = [];

    fprintf('%f length  res_geoid_err_w\n',length (Grid_res_geoid_err_w));
    fprintf('%f min     res_geoid_err_w\n',min    (Grid_res_geoid_err_w));
    fprintf('%f max     res_geoid_err_w\n',max    (Grid_res_geoid_err_w));
    fprintf('%f mean    res_geoid_err_w\n',mean   (Grid_res_geoid_err_w));
    fprintf('%f median  res_geoid_err_w\n',median (Grid_res_geoid_err_w));
    fprintf('%f std     res_geoid_err_w\n',std    (Grid_res_geoid_err_w));

    
%     fprintf('%f length  GridResGravErrW\n',length (GridResGravErrW));
%     fprintf('%f min     GridResGravErrW\n',min    (GridResGravErrW));
%     fprintf('%f max     GridResGravErrW\n',max    (GridResGravErrW));
%     fprintf('%f mean    GridResGravErrW\n',mean   (GridResGravErrW));
%     fprintf('%f median  GridResGravErrW\n',median (GridResGravErrW));
%     fprintf('%f std     GridResGravErrW\n',std    (GridResGravErrW));
    
end
