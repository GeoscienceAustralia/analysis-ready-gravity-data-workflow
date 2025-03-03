function plotMosaicTiles(Coastline,GRID_PARA,LongDEM,LatDEM,Grid_res_geoid_w,resAGQG,ZDeg,Lev,Vals_Lev, AGQG_Vals_Lev, ...
    Grid_res_geoid_err_w,Grid_res_grav_w,Grid_res_grav_Bouguer_w,Grid_res_grav_err_w,plotsFolder)

    % common variables for plotting
    axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
    axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
    axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
    axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
    axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;

    % plot residualGeoidvsAGQG
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    subplot(1,2,1);
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid_w)
    customizeMap('residualGeoid','m',Coastline,axisLimits)
    caxis([-0.5 0.5])
     
    subplot(1,2,2);
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),resAGQG-ZDeg)
    customizeMap('residualAGQG','m',Coastline,axisLimits) 
    caxis([-0.5 0.5])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualGeoidvsAGQG','.png']) 

    % plot GPSlevelling vs LSC
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    scatter(Lev(:,1),Lev(:,2),15,Vals_Lev-mean(Vals_Lev(~isnan(Vals_Lev))),'filled')
    customizeMap('Geometric and LSC Gravimetric Geoid Difference','m',Coastline,axisLimits)
    caxis([-0.25 0.25])
    %saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevellingLSC','.png']) 
    
    % plot GPSlevelling vs reference AGQG
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    scatter(Lev(:,1),Lev(:,2),15,AGQG_Vals_Lev-mean(AGQG_Vals_Lev(~isnan(AGQG_Vals_Lev))),'filled')
    customizeMap('Geometric and AGQG Difference','m',Coastline,axisLimits)
    caxis([-0.25 0.25])
    %saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevellingAGQG','.png']) 

    % plot difference LSC and AGQG at GPSlevelling
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    scatter(Lev(:,1),Lev(:,2),15,AGQG_Vals_Lev-mean(AGQG_Vals_Lev(~isnan(AGQG_Vals_Lev)))- ...
        Vals_Lev+mean(Vals_Lev(~isnan(Vals_Lev))),'filled')
    customizeMap('AGQG and LSC Gravimetric Geoid Difference','m',Coastline,axisLimits)
    caxis([-0.25 0.25])
    %saveas(gcf,[plotsFolder,'MosaicTiles','AGQGvsLSCdiff','.png']) 
 
    % plot residualGeoidWeighted-residualAGQG
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid_w-(resAGQG-ZDeg))
    customizeMap('ResidualDiffGeoidAndAGQG','m',Coastline,axisLimits)
    caxis([-1 1])
    saveas(gcf,[plotsFolder,'MosaicTiles','ResidualDiffGeoidAndAGQG','.png'])
 
    % plot residualGeoidSigmaError
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid_err_w)
    customizeMap('residualGeoidSigmaError','m',Coastline,axisLimits)
    %caxis([0 0.075])%
    %caxis([0 0.03])
    caxis([0 0.09])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualGeoidSigmaError','.png'])
    
    % plot residualGravityWeighted
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_grav_w)
    customizeMap('residualFreeAirGravityWeighted','mGal',Coastline,axisLimits)
    saveas(gcf,[plotsFolder,'MosaicTiles','residualFreeAirGravityWeighted','.png'])
 
    % plot residualBouguerGravityWeighted
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_grav_Bouguer_w)
    customizeMap('residualBouguerGravityWeighted','mGal',Coastline,axisLimits)
    %caxis([-110 110])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualBouguerGravityWeighted','.png'])
 
    % plot residualGravitySigmaError
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_grav_err_w)
    customizeMap('residualGravitySigmaError','mGal',Coastline,axisLimits)
    %caxis([0 55])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualGravitySigmaError','.png'])
end
