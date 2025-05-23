function plotOutputData(LongDEM,LatDEM,Grid_res_geoid,Weights,Gravity,Gravity_BA,Grid_res_grav,GGM_Zeta_griddedInterpolant, ...
    Grid_res_geoid_err,Grid_res_grav_err,Lev,Coastline,GRID_PARA,LEVELLING_PARA,block_counter,LONGsi,LATsi)

        % common variables for plotting
        axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
        axisLimits.lonMinLimit=LONGsi-GRID_PARA.buffer;
        axisLimits.lonMaxLimit=LONGsi+GRID_PARA.buffer;
        axisLimits.latMinLimit=LATsi-GRID_PARA.buffer;
        axisLimits.latMaxLimit=LATsi+GRID_PARA.buffer;

        % plot Residual Geoid 
        figure('Name','computeLSC','NumberTitle','off'); 
        clf
        hold on
        imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid./Weights)
        customizeMap('m',Coastline,axisLimits)
        title(['Residual Geoid for Block ',num2str(block_counter),' Centred at: ',num2str(LONGsi),', ',num2str(LATsi)])
        caxis([-2*std(Grid_res_geoid(:)./Weights(:),"omitnan") 2*std(Grid_res_geoid(:)./Weights(:),"omitnan")])
        drawnow
        %saveas(gcf,['Outputs/plots/ResidualGeoid',num2str(block_counter),'.png'])

        % plot Free-air Gravity Anomaly
        figure('Name','computeLSC','NumberTitle','off'); 
        clf
        hold on
        imagesc(LongDEM(1,:),LatDEM(:,1),Gravity)
        customizeMap('mGal',Coastline,axisLimits)
        title(['Free-air Gravity Anomaly for Block ',num2str(block_counter),' Centred at: ',num2str(LONGsi),', ',num2str(LATsi)])
        caxis([-2*std(Grid_res_grav(:)./Weights(:),"omitnan") 2*std(Grid_res_grav(:)./Weights(:),"omitnan")])
        drawnow
        %saveas(gcf,['Outputs/plots/FreeAirGravity',num2str(block_counter),'.png'])

        % plot LSC Geoid
        figure('Name','computeLSC','NumberTitle','off'); 
        clf
        hold on
        imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid./Weights+GGM_Zeta_griddedInterpolant(LongDEM,-LatDEM,LatDEM*0))
        customizeMap('m',Coastline,axisLimits)
        title(['LSC Geoid for Block ',num2str(block_counter),' Centred at: ',num2str(LONGsi),', ',num2str(LATsi)])
        drawnow
        %saveas(gcf,['Outputs/plots/LSCGeoid',num2str(block_counter),'.png'])

        % plot LSC Geoid Error
        figure('Name','computeLSC','NumberTitle','off'); 
        clf
        hold on
        imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid_err./Weights)
        customizeMap('m',Coastline,axisLimits)
        if LEVELLING_PARA.Plot_Stats && LEVELLING_PARA.Lev_eval    
        plot(Lev(:,1),Lev(:,2),'g*')
        end
        caxis([0 min(max(Grid_res_geoid_err(:)./Weights(:)),0.1)])
        title(['LSC Geoid Error for Block ',num2str(block_counter),' Centred at: ',num2str(LONGsi),', ',num2str(LATsi)])
        drawnow
        %saveas(gcf,['Outputs/plots/LSCGeoidError',num2str(block_counter),'.png'])

        % plot Gravity Error
        figure('Name','computeLSC','NumberTitle','off'); 
        clf
        hold on
        imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_grav_err./Weights)
        customizeMap('mGal',Coastline,axisLimits)
        caxis([0 min(max(Grid_res_grav_err(:)./Weights(:)),10)])
        title(['Gravity Error for Block ',num2str(block_counter),' Centred at: ',num2str(LONGsi),', ',num2str(LATsi)])
        drawnow
        %saveas(gcf,['Outputs/plots/GravityError',num2str(block_counter),'.png'])

        % plot Bouguer Gravity Anomaly
        figure('Name','computeLSC','NumberTitle','off'); 
        clf
        hold on
        imagesc(LongDEM(1,:),LatDEM(:,1),Gravity_BA)
        customizeMap('mGal',Coastline,axisLimits)
        title(['Bouguer Gravity Anomaly for Block ',num2str(block_counter),' Centred at: ',num2str(LONGsi),', ',num2str(LATsi)])
        drawnow
        %saveas(gcf,['Outputs/plots/BouguerGravityAnomaly',num2str(block_counter),'.png'])
end



       
        