function plotCovParameters(covParameters,GRID_PARA,Coastline)

        % common variables for plotting
        axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
        axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
        axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
        axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
        axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;
        
        % plots
        figure('Name','MosaicTiles','NumberTitle','off'); 
        hold on
        scatter (covParameters.filenameLong,-covParameters.filenameLat,15,covParameters.COV_PARA_RTM_A,'filled')
        customizeMap('COV PARA RTM A','',Coastline,axisLimits)
        %caxis([0 220])
        
        figure('Name','MosaicTiles','NumberTitle','off'); 
        hold on
        scatter (covParameters.filenameLong,-covParameters.filenameLat,15,covParameters.COV_PARA_RTM_B,'filled')
        customizeMap('COV PARA RTM B','',Coastline,axisLimits)
        
        figure('Name','MosaicTiles','NumberTitle','off'); 
        hold on
        scatter (covParameters.filenameLong,-covParameters.filenameLat,15,covParameters.COV_PARA_RTM_M,'filled')
        customizeMap('COV PARA RTM M','',Coastline,axisLimits)
       
        figure('Name','MosaicTiles','NumberTitle','off'); 
        hold on
        scatter (covParameters.filenameLong,-covParameters.filenameLat,15,covParameters.COV_PARA_Faye_A,'filled')
        customizeMap('COV PARA Faye A','',Coastline,axisLimits)
        %caxis([0 220])
        
        figure('Name','MosaicTiles','NumberTitle','off'); 
        hold on
        scatter (covParameters.filenameLong,-covParameters.filenameLat,15,covParameters.COV_PARA_Faye_B,'filled')
        customizeMap('COV PARA Faye B','',Coastline,axisLimits)
        
        figure('Name','MosaicTiles','NumberTitle','off'); 
        hold on
        scatter (covParameters.filenameLong,-covParameters.filenameLat,15,covParameters.COV_PARA_Faye_M,'filled')
        customizeMap('COV PARA Faye M','',Coastline,axisLimits)

end




 

