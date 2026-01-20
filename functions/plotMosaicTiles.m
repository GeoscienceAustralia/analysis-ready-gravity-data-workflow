function plotMosaicTiles(Coastline,GRID_PARA,LongDEM,LatDEM,Grid_res_geoid_w,resAGQG,ZDeg,Lev,geoGravGeoidDiff, geoRefGravGeoidDiff, ...
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
    customizeMap('Residual LSC AGQG','m',Coastline,axisLimits)
    caxis([-0.5 0.5])
     
    subplot(1,2,2);
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),resAGQG-ZDeg)
    customizeMap('Residual 2022 AGQG','m',Coastline,axisLimits) 
    caxis([-0.5 0.5])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualLSCvs2022AGQG','.png']) 

    % plot GPSlevelling vs LSC
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    scatter(Lev(:,1),Lev(:,2),15,geoGravGeoidDiff-mean(geoGravGeoidDiff(~isnan(geoGravGeoidDiff))),'filled')
    customizeMap('Geometric and LSC AGQG Difference','m',Coastline,axisLimits)
    caxis([-0.3 0.3])
    saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevellingLSCAGQG','.fig']) 
    saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevellingLSCAGQG','.png']) 

    validVals = geoGravGeoidDiff(~isnan(geoGravGeoidDiff));
    valDiff = geoGravGeoidDiff - mean(validVals);
    
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
    scatter(Lev(:,1),Lev(:,2),15,geoRefGravGeoidDiff-mean(geoRefGravGeoidDiff(~isnan(geoRefGravGeoidDiff))),'filled')
    customizeMap('Geometric and 2022 AGQG Difference','m',Coastline,axisLimits)
    caxis([-0.3 0.3])
    saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevelling2022AGQG','.fig']) 
    saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevelling2022AGQG','.png'])

    validVals = geoRefGravGeoidDiff(~isnan(geoRefGravGeoidDiff));
    valDiff = geoRefGravGeoidDiff - mean(validVals);
    
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
    scatter(Lev(:,1),Lev(:,2),15,geoRefGravGeoidDiff-mean(geoRefGravGeoidDiff(~isnan(geoRefGravGeoidDiff)))- ...
        geoGravGeoidDiff+mean(geoGravGeoidDiff(~isnan(geoGravGeoidDiff))),'filled')
    customizeMap('2022 AGQG and LSC AGQG Difference','m',Coastline,axisLimits)
    caxis([-0.1 0.1])
    saveas(gcf,[plotsFolder,'MosaicTiles','oldVSnewAGQGdiffGPSpoints','.fig']) 
    saveas(gcf,[plotsFolder,'MosaicTiles','oldVSnewAGQGdiffGPSpoints','.png']) 
 
    % plot residualGeoidWeighted-residualAGQG
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),(resAGQG-ZDeg-Grid_res_geoid_w))
    customizeMap('2022 AGQG and LSC AGQG Difference','m',Coastline,axisLimits)
    caxis([-0.1 0.1])
    saveas(gcf,[plotsFolder,'MosaicTiles','DiffResidualLSC2022AGQG','.png'])
 
    % plot residualGeoidSigmaError
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    imagesc(LongDEM(1,:),LatDEM(:,1),Grid_res_geoid_err_w)
    customizeMap('Residual AGQG Sigma Error','m',Coastline,axisLimits)
    caxis([0 0.03])%caxis([0 0.09])
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
    caxis([0 10])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualGravitySigmaError','.png'])
    saveas(gcf,[plotsFolder,'MosaicTiles','residualGravitySigmaError','.fig'])
end
