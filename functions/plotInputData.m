function plotInputData(Gravo,Coastline,GRID_PARA)
    % Plot input data: 
    % be carefull, the size of AusDEM is 15462841
    %plotCustomScatter(DEM_data(:,1),DEM_data(:,2),DEM_data(:,3),GRID_PARA,'DEM','m',Coastline,[],OUTPUT_PARA.plotsFolder)
    
    plotCustomScatter(Gravo(:,1),Gravo(:,2),Gravo(:,3),GRID_PARA,'GravityTopographyHeight','m',Coastline,[0 2000],OUTPUT_PARA.plotsFolder)
    
    plotCustomScatter(Gravo(:,1),Gravo(:,2),Gravo(:,4),GRID_PARA,'Gravity','mGal',Coastline,[],OUTPUT_PARA.plotsFolder)
    
    plotCustomScatter(Gravo(:,1),Gravo(:,2),Gravo(:,5),GRID_PARA,'GravityUncertainty','mGal',Coastline,[],OUTPUT_PARA.plotsFolder)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                                                                                    
    %plotCustomScatter(gravGradFiltered(:,1),gravGradFiltered(:,2),gravGradFiltered(:,4),GRID_PARA,'GravityGradient','mGal/m',Coastline,[],OUTPUT_PARA.plotsFolder)
    
    %plotProfiles(gravGradFiltered(1:180,2),gravGradFiltered(1:180,3),gravGradFiltered(1:180,4),gravGradFiltered(1,1),'Latitude','GravityGradientFlightAltitude [m]','GravityGradient [mGal/m]','GravityGradientProfile1',OUTPUT_PARA.plotsFolder)
    
    %plotProfiles(gravGradFiltered(180+1:2*180,2),gravGradFiltered(180+1:2*180,3),gravGradFiltered(180+1:2*180,4),gravGradFiltered(180+1,1),'Latitude','GravityGradientFlightAltitude [m]','GravityGradient [mGal/m]','GravityGradientProfile2',OUTPUT_PARA.plotsFolder)
end