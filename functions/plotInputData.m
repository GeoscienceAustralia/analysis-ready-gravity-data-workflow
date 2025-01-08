function plotInputData(Gravo,Coastline,GRID_PARA,OUTPUT_PARA)
    % Plot input data: 
    % be carefull, the size of AusDEM is 15462841
    % plotCustomScatter(DEM_data(:,1),DEM_data(:,2),DEM_data(:,3),GRID_PARA,'DEM','m',Coastline,[],OUTPUT_PARA.plotsFolder)
    % for perth simulation 
    % Define the coordinates of the polygon
    lon = [115.4333, 116.0500, 116.2500, 116.2500, 115.6167, 115.6167];
    lat = [-31.4500, -31.4500, -32.0000, -32.5833, -32.5833, -32.0000];
    
    plotCustomScatter(Gravo(:,1),Gravo(:,2),Gravo(:,3),GRID_PARA,'GravityTopographyHeight','m',Coastline,[],OUTPUT_PARA.plotsFolder)%[0 2000]
    hold on;
    plot(lon, lat, 'b-', 'LineWidth', 2); % Plot the polygon
    hold on;
    plot([lon lon(1)], [lat lat(1)], 'b-', 'LineWidth', 2); % Close the polygon

    plotCustomScatter(Gravo(:,1),Gravo(:,2),Gravo(:,4),GRID_PARA,'Gravity','mGal',Coastline,[],OUTPUT_PARA.plotsFolder)
    hold on;
    plot(lon, lat, 'b-', 'LineWidth', 2); % Plot the polygon
    hold on;
    plot([lon lon(1)], [lat lat(1)], 'b-', 'LineWidth', 2); % Close the polygon
   
    plotCustomScatter(Gravo(:,1),Gravo(:,2),Gravo(:,5),GRID_PARA,'GravityUncertainty','mGal',Coastline,[],OUTPUT_PARA.plotsFolder)
    hold on;
    plot(lon, lat, 'b-', 'LineWidth', 2); % Plot the polygon
    hold on;
    plot([lon lon(1)], [lat lat(1)], 'b-', 'LineWidth', 2); % Close the polygon
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                                                                                    
    %plotCustomScatter(gravGradFiltered(:,1),gravGradFiltered(:,2),gravGradFiltered(:,4),GRID_PARA,'GravityGradient','mGal/m',Coastline,[],OUTPUT_PARA.plotsFolder)
    
    %plotProfiles(gravGradFiltered(1:180,2),gravGradFiltered(1:180,3),gravGradFiltered(1:180,4),gravGradFiltered(1,1),'Latitude','GravityGradientFlightAltitude [m]','GravityGradient [mGal/m]','GravityGradientProfile1',OUTPUT_PARA.plotsFolder)
    
    %plotProfiles(gravGradFiltered(180+1:2*180,2),gravGradFiltered(180+1:2*180,3),gravGradFiltered(180+1:2*180,4),gravGradFiltered(180+1,1),'Latitude','GravityGradientFlightAltitude [m]','GravityGradient [mGal/m]','GravityGradientProfile2',OUTPUT_PARA.plotsFolder)
end