function plotInputData(Gravo, gravGradFiltered, Coastline, GRID_PARA, OUTPUT_PARA,DEMdata)
    % Plot input data: 
    % be careful, the size of AusDEM is 15462841
    %plotCustomScatter(DEMdata(:,1),DEMdata(:,2),DEMdata(:,3),GRID_PARA,'DEM','m',Coastline,[],OUTPUT_PARA.plotsFolder)
    % for Perth simulation 

    % Plot GravityTopographyHeight
    plotCustomScatter(Gravo(:,1), Gravo(:,2), Gravo(:,3), GRID_PARA, 'Height', 'm', Coastline, [], OUTPUT_PARA.plotsFolder);
    % If lon and lat are not empty, plot the polygon
    if ~isempty(OUTPUT_PARA.polygonLon) && ~isempty(OUTPUT_PARA.polygonLat)
        hold on;
        plot(OUTPUT_PARA.polygonLon, OUTPUT_PARA.polygonLat, 'magenta-', 'LineWidth', 2); % Plot the polygon
        hold off; % Release the hold to save the plot with the polygon
        % Save the figure
        saveas(gcf,[OUTPUT_PARA.plotsFolder,'scatter','Height','.png'])
    end

    % Plot Gravity
    plotCustomScatter(Gravo(:,1), Gravo(:,2), Gravo(:,4), GRID_PARA, 'Gravity', 'mGal', Coastline, [], OUTPUT_PARA.plotsFolder);
    % If lon and lat are not empty, plot the polygon
    if ~isempty(OUTPUT_PARA.polygonLon) && ~isempty(OUTPUT_PARA.polygonLat)
        hold on;
        plot(OUTPUT_PARA.polygonLon, OUTPUT_PARA.polygonLat, 'magenta-', 'LineWidth', 2); % Plot the polygon
        hold off; % Release the hold to save the plot with the polygon
        % Save the figure
        saveas(gcf,[OUTPUT_PARA.plotsFolder,'scatter','Gravity','.png'])
    end

    % Plot GravityUncertainty
    plotCustomScatter(Gravo(:,1), Gravo(:,2), Gravo(:,5), GRID_PARA, 'Uncertainty', 'mGal', Coastline, [], OUTPUT_PARA.plotsFolder);
    % If lon and lat are not empty, plot the polygon
    if ~isempty(OUTPUT_PARA.polygonLon) && ~isempty(OUTPUT_PARA.polygonLat)
        hold on;
        plot(OUTPUT_PARA.polygonLon, OUTPUT_PARA.polygonLat, 'magenta-', 'LineWidth', 2); % Plot the polygon
        hold off; % Release the hold to save the plot with the polygon
        % Save the figure
        saveas(gcf,[OUTPUT_PARA.plotsFolder,'scatter','Uncertainty','.png'])
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(gravGradFiltered)
                                                                                                                                                                    
    plotCustomScatter(gravGradFiltered(:,1),gravGradFiltered(:,2),gravGradFiltered(:,4),GRID_PARA,'Gravity Gradient','mGal/m',Coastline,[],OUTPUT_PARA.plotsFolder)
    plotCustomScatter(gravGradFiltered(:,1),gravGradFiltered(:,2),gravGradFiltered(:,5),GRID_PARA,'Gravity Gradient Uncertainty','mGal/m',Coastline,[],OUTPUT_PARA.plotsFolder)
    plotCustomScatter(gravGradFiltered(:,1),gravGradFiltered(:,2),gravGradFiltered(:,3),GRID_PARA,'Gravity Gradient Height','m',Coastline,[],OUTPUT_PARA.plotsFolder)
    %plotProfiles(gravGradFiltered(1:180,2),gravGradFiltered(1:180,3),gravGradFiltered(1:180,4),gravGradFiltered(1,1),'Latitude','GravityGradientFlightAltitude [m]','GravityGradient [mGal/m]','GravityGradientProfile1',OUTPUT_PARA.plotsFolder)
    %plotProfiles(gravGradFiltered(180+1:2*180,2),gravGradFiltered(180+1:2*180,3),gravGradFiltered(180+1:2*180,4),gravGradFiltered(180+1,1),'Latitude','GravityGradientFlightAltitude [m]','GravityGradient [mGal/m]','GravityGradientProfile2',OUTPUT_PARA.plotsFolder)
    end
end