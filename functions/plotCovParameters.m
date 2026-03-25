function plotCovParameters(covParameters,GRID_PARA,Coastline,plotsFolder)

        % common variables for plotting
        axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
        axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
        axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
        axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
        axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;

        % disply statistics
        fprintf('%-8s %12s %12s %12s %12s %12s\n', ...
            'Set','Min','Max','Median','Mean','Std');
        
        printStats = @(name,x) fprintf( ...
            '%-8s %12.6g %12.6g %12.6g %12.6g %12.6g\n', ...
            name, ...
            min(x,[], 'omitnan'), ...
            max(x,[], 'omitnan'), ...
            median(x,  'omitnan'), ...
            mean(x,    'omitnan'), ...
            std(x, 0,  'omitnan'));
        
        printStats('RTM_A',  covParameters.COV_PARA_RTM_A)
        printStats('RTM_B',  covParameters.COV_PARA_RTM_B)
        printStats('Faye_A', covParameters.COV_PARA_Faye_A)
        printStats('Faye_B', covParameters.COV_PARA_Faye_B)

        % 1D plots
        figure('Name','MosaicTiles','NumberTitle','off'); 
        subplot(3,1,1)
        plot(covParameters.COV_PARA_RTM_A,'b.')
        title('COVPARA RTM A')
        subplot(3,1,2)
        plot(covParameters.COV_PARA_RTM_B,'g.')
        title('COVPARA RTM B')
        subplot(3,1,3)
        plot(covParameters.COV_PARA_RTM_M,'m.')
        title('COVPARA RTM M')
        saveas(gcf,[plotsFolder,'MosaicTiles','COVPARARTM','.png']) 
        
        figure('Name','MosaicTiles','NumberTitle','off'); 
        subplot(3,1,1)
        plot(covParameters.COV_PARA_Faye_A,'b.')
        title('COVPARA Faye A')
        subplot(3,1,2)
        plot(covParameters.COV_PARA_Faye_B,'g.')
        title('COVPARA Faye B')
        subplot(3,1,3)
        plot(covParameters.COV_PARA_Faye_M,'m.')
        title('COVPAR Faye M')
        saveas(gcf,[plotsFolder,'MosaicTiles','COVPARARFaye','.png']) 

        % 2D plots
        figure('Name','MosaicTiles','NumberTitle','off'); 
        hold on
        scatter (covParameters.filenameLong,-covParameters.filenameLat,15,covParameters.COV_PARA_RTM_A,'filled')
        customizeMap('COV PARA RTM A','',Coastline,axisLimits)
        caxis([1 30])
        saveas(gcf,[plotsFolder,'MosaicTiles','COVPARARTMA','.png']) 
        
        figure('Name','MosaicTiles','NumberTitle','off'); 
        hold on
        scatter (covParameters.filenameLong,-covParameters.filenameLat,15,covParameters.COV_PARA_RTM_B,'filled')
        customizeMap('COV PARA RTM B','',Coastline,axisLimits)
        caxis([-300 -200])
        saveas(gcf,[plotsFolder,'MosaicTiles','COVPARARTMB','.png'])
        
        figure('Name','MosaicTiles','NumberTitle','off'); 
        hold on
        scatter (covParameters.filenameLong,-covParameters.filenameLat,15,covParameters.COV_PARA_RTM_M,'filled')
        customizeMap('COV PARA RTM M','',Coastline,axisLimits)
        saveas(gcf,[plotsFolder,'MosaicTiles','COVPARARTMM','.png'])
       
        figure('Name','MosaicTiles','NumberTitle','off'); 
        hold on
        scatter (covParameters.filenameLong,-covParameters.filenameLat,15,covParameters.COV_PARA_Faye_A,'filled')
        customizeMap('COV PARA Faye A','',Coastline,axisLimits)
        caxis([1 30])
        saveas(gcf,[plotsFolder,'MosaicTiles','COVPARAFayeA','.png'])

        figure('Name','MosaicTiles','NumberTitle','off'); 
        hold on
        scatter (covParameters.filenameLong,-covParameters.filenameLat,15,covParameters.COV_PARA_Faye_B,'filled')
        customizeMap('COV PARA Faye B','',Coastline,axisLimits)
       caxis([-300 -200])
        saveas(gcf,[plotsFolder,'MosaicTiles','COVPARAFayeB','.png'])

        figure('Name','MosaicTiles','NumberTitle','off'); 
        hold on
        scatter (covParameters.filenameLong,-covParameters.filenameLat,15,covParameters.COV_PARA_Faye_M,'filled')
        customizeMap('COV PARA Faye M','',Coastline,axisLimits)
        saveas(gcf,[plotsFolder,'MosaicTiles','COVPARAFayeM','.png'])
end




 

