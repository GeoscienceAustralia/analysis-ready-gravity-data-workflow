function DisplayAreaStatistics(Coastline,GRID_PARA,LongDEM,LatDEM,Grid_res_geoid_w, ...
    Grid_res_geoid_err_w,OUTPUT_PARA)

    % common variables for plotting
    axisLimits.latMeanCosine=abs(cos(deg2rad(mean([min(OUTPUT_PARA.polygonLat) max(OUTPUT_PARA.polygonLat)]))));
    axisLimits.lonMinLimit=min(OUTPUT_PARA.polygonLon)-GRID_PARA.buffer;
    axisLimits.lonMaxLimit=max(OUTPUT_PARA.polygonLon)+GRID_PARA.buffer;
    axisLimits.latMinLimit=min(OUTPUT_PARA.polygonLat)-GRID_PARA.buffer;
    axisLimits.latMaxLimit=max(OUTPUT_PARA.polygonLat)+GRID_PARA.buffer;

    if ~isempty(OUTPUT_PARA.polygonLon) && ~isempty(OUTPUT_PARA.polygonLat)

       areaIn=inpolygon(LongDEM,LatDEM,OUTPUT_PARA.polygonLon,OUTPUT_PARA.polygonLat);
       Grid_res_geoid_err_w(areaIn==0)=nan;
       Grid_res_geoid_w(areaIn==0)=nan;
       nanMask = ~isnan(Grid_res_geoid_err_w);% true where data are valid

       if OUTPUT_PARA.PLOT_GRIDS

            figure('Name','MosaicTiles','NumberTitle','off'); 
            clf
            subplot(1,2,1);
            hold on
            imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid_w,'AlphaData', nanMask)
            customizeMap('Residual LSC AGQG','m',Coastline,axisLimits)
            caxis([-0.5 0.5])
           
            subplot(1,2,2);
            hold on
            imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid_err_w,'AlphaData', nanMask)
            %plot(OUTPUT_PARA.polygonLon, OUTPUT_PARA.polygonLat, 'magenta-', 'LineWidth', 2); % Plot the polygon
            customizeMap('Residual AGQG Sigma Error','m',Coastline,axisLimits)
            caxis([0 0.02])
            saveas(gcf,[OUTPUT_PARA.plotsFolder,'MosaicTiles','residualGeoidFocused','.png'])

        end

    end

    disp('Display statistics')

    Grid_res_geoid_err_w(isnan(Grid_res_geoid_err_w)) = [];

    fprintf('%f length  res_geoid_err_w\n',length (Grid_res_geoid_err_w));
    fprintf('%f min     res_geoid_err_w\n',min    (Grid_res_geoid_err_w));
    fprintf('%f max     res_geoid_err_w\n',max    (Grid_res_geoid_err_w));
    fprintf('%f mean    res_geoid_err_w\n',mean   (Grid_res_geoid_err_w));
    fprintf('%f median  res_geoid_err_w\n',median (Grid_res_geoid_err_w));
    fprintf('%f std     res_geoid_err_w\n',std    (Grid_res_geoid_err_w));

    Grid_res_geoid_w(isnan(Grid_res_geoid_w)) = [];

    fprintf('%f length  res_geoid_w\n',length (Grid_res_geoid_w));
    fprintf('%f min     res_geoid_w\n',min    (Grid_res_geoid_w));
    fprintf('%f max     res_geoid_w\n',max    (Grid_res_geoid_w));
    fprintf('%f mean    res_geoid_w\n',mean   (Grid_res_geoid_w));
    fprintf('%f median  res_geoid_w\n',median (Grid_res_geoid_w));
    fprintf('%f std     res_geoid_w\n',std    (Grid_res_geoid_w));

    
%     fprintf('%f length  GridResGravErrW\n',length (GridResGravErrW));
%     fprintf('%f min     GridResGravErrW\n',min    (GridResGravErrW));
%     fprintf('%f max     GridResGravErrW\n',max    (GridResGravErrW));
%     fprintf('%f mean    GridResGravErrW\n',mean   (GridResGravErrW));
%     fprintf('%f median  GridResGravErrW\n',median (GridResGravErrW));
%     fprintf('%f std     GridResGravErrW\n',std    (GridResGravErrW));
    
end
