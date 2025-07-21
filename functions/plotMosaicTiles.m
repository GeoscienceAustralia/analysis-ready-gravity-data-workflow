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
    customizeMap('Residual Geoid','m',Coastline,axisLimits)
    caxis([-0.5 0.5])
     
    subplot(1,2,2);
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),resAGQG-ZDeg)
    customizeMap('Residual AGQG','m',Coastline,axisLimits) 
    caxis([-0.5 0.5])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualGeoidvsAGQG','.png']) 

    % plot GPSlevelling vs LSC
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    scatter(Lev(:,1),Lev(:,2),15,Vals_Lev-mean(Vals_Lev(~isnan(Vals_Lev))),'filled')
    customizeMap('Geometric and LSC Gravimetric Geoid Difference','m',Coastline,axisLimits)
    caxis([-0.1 0.1])
    saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevellingLSC','.fig']) 

    validVals = Vals_Lev(~isnan(Vals_Lev));
    valDiff = Vals_Lev - mean(validVals);
    
    fprintf('%f length  GPSlevellingLSC\n',length (valDiff));
    fprintf('%f min     GPSlevellingLSC\n',min    (valDiff));
    fprintf('%f max     GPSlevellingLSC\n',max    (valDiff));
    fprintf('%f mean    GPSlevellingLSC\n',mean   (valDiff));
    fprintf('%f median  GPSlevellingLSC\n',median (valDiff));
    fprintf('%f std     GPSlevellingLSC\n',std    (valDiff));

    % plot GPSlevelling vs reference AGQG
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    scatter(Lev(:,1),Lev(:,2),15,AGQG_Vals_Lev-mean(AGQG_Vals_Lev(~isnan(AGQG_Vals_Lev))),'filled')
    customizeMap('Geometric and AGQG Difference','m',Coastline,axisLimits)
    caxis([-0.1 0.1])
    saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevellingAGQG','.fig']) 

    validVals = AGQG_Vals_Lev(~isnan(AGQG_Vals_Lev));
    valDiff = AGQG_Vals_Lev - mean(validVals);
    
    fprintf('%f length  GPSlevellingAGQG\n',length (valDiff));
    fprintf('%f min     GPSlevellingAGQG\n',min    (valDiff));
    fprintf('%f max     GPSlevellingAGQG\n',max    (valDiff));
    fprintf('%f mean    GPSlevellingAGQG\n',mean   (valDiff));
    fprintf('%f median  GPSlevellingAGQG\n',median (valDiff));
    fprintf('%f std     GPSlevellingAGQG\n',std    (valDiff));

    % plot difference LSC and AGQG at GPSlevelling
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    scatter(Lev(:,1),Lev(:,2),15,AGQG_Vals_Lev-mean(AGQG_Vals_Lev(~isnan(AGQG_Vals_Lev)))- ...
        Vals_Lev+mean(Vals_Lev(~isnan(Vals_Lev))),'filled')
    customizeMap('AGQG and LSC Gravimetric Geoid Difference','m',Coastline,axisLimits)
    caxis([-0.1 0.1])
    saveas(gcf,[plotsFolder,'MosaicTiles','AGQGvsLSCdiff','.fig']) 
 
    % plot residualGeoidWeighted-residualAGQG
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid_w-(resAGQG-ZDeg))
    customizeMap('Residual Diff Geoid and AGQG','m',Coastline,axisLimits)
    caxis([-1 1])
    saveas(gcf,[plotsFolder,'MosaicTiles','ResidualDiffGeoidAndAGQG','.png'])
 
    % plot residualGeoidSigmaError
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid_err_w)
    customizeMap('Residual Geoid Sigma Error','m',Coastline,axisLimits)
    %caxis([0 0.075])%caxis([0 0.03])%caxis([0 0.09])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualGeoidSigmaError','.png'])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualGeoidSigmaError','.fig'])
    
    % plot residualGravityWeighted
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_grav_w)
    customizeMap('Residual Free Air Gravity Weighted','mGal',Coastline,axisLimits)
    saveas(gcf,[plotsFolder,'MosaicTiles','residualFreeAirGravityWeighted','.png'])
 
    % plot residualBouguerGravityWeighted
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_grav_Bouguer_w)
    customizeMap('Residual Bouguer Gravity Weighted','mGal',Coastline,axisLimits)
    %caxis([-110 110])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualBouguerGravityWeighted','.png'])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualBouguerGravityWeighted','.fig'])
 
    % plot residualGravitySigmaError
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_grav_err_w)
    customizeMap('Residual Gravity Sigma Error','mGal',Coastline,axisLimits)
    %caxis([0 55])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualGravitySigmaError','.png'])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualGravitySigmaError','.fig'])
end
