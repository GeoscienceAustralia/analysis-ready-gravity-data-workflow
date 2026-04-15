function plotGPSlevelling(Coastline,GRID_PARA,Lev,geoGravGeoidDiff, geoRefGravGeoidDiff,plotsFolder)

    % common variables for plottinggeoGravGeoidDiff
    axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
    axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
    axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
    axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
    axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;
    
    % plot GPSlevelling vs LSC
    validVals = geoGravGeoidDiff(~isnan(geoGravGeoidDiff));
    valDiff = geoGravGeoidDiff - mean(validVals);
    
    figure('Name','geometric','NumberTitle','off'); 
    clf
    hold on
    scatter(Lev(:,1),Lev(:,2),5,valDiff,'filled')
    customizeMap('Geometric and LSC AGQG Difference','m',Coastline,axisLimits)
    caxis([-0.3 0.3])
    saveas(gcf,[plotsFolder,'geometric','GPSlevellingLSCAGQG','.fig']) 
    saveas(gcf,[plotsFolder,'geometric','GPSlevellingLSCAGQG','.png']) 
    
    fprintf('Geometric and LSC AGQG Difference Statistics:\n');
    fprintf('  Count: %d\n',     numel(valDiff));
    fprintf('  Mean: %.4f\n',    mean(valDiff));
    fprintf('  Median: %.4f\n',  median(valDiff));
    fprintf('  Std Dev: %.4f\n', std(valDiff));
    fprintf('  Min: %.4f\n',     min(valDiff));
    fprintf('  Max: %.4f\n\n',   max(valDiff));
    
    % plot GPSlevelling vs reference AGQG
    
    validValsRef = geoRefGravGeoidDiff(~isnan(geoRefGravGeoidDiff));
    valDiffRef = geoRefGravGeoidDiff - mean(validValsRef);
    
    figure('Name','geometric','NumberTitle','off'); 
    clf
    hold on
    scatter(Lev(:,1),Lev(:,2),5,valDiffRef,'filled')
    customizeMap('Geometric and 2022 AGQG Difference','m',Coastline,axisLimits)
    caxis([-0.3 0.3])
    saveas(gcf,[plotsFolder,'geometric','GPSlevelling2022AGQG','.fig']) 
    saveas(gcf,[plotsFolder,'geometric','GPSlevelling2022AGQG','.png'])
    
    fprintf('Geometric and 2022 AGQG Difference Statistics:\n');
    fprintf('  Count: %d\n',     numel(valDiffRef));
    fprintf('  Mean: %.4f\n',    mean(valDiffRef));
    fprintf('  Median: %.4f\n',  median(valDiffRef));
    fprintf('  Std Dev: %.4f\n', std(valDiffRef));
    fprintf('  Min: %.4f\n',     min(valDiffRef));
    fprintf('  Max: %.4f\n\n',   max(valDiffRef));
    
    % plot difference LSC and AGQG at GPSlevelling
    figure('Name','geometric','NumberTitle','off'); 
    clf
    hold on
    scatter(Lev(:,1),Lev(:,2),5,geoRefGravGeoidDiff-mean(geoRefGravGeoidDiff(~isnan(geoRefGravGeoidDiff)))- ...
        geoGravGeoidDiff+mean(geoGravGeoidDiff(~isnan(geoGravGeoidDiff))),'filled')
    customizeMap('2022 AGQG and LSC AGQG Difference','m',Coastline,axisLimits)
    caxis([-0.1 0.1])
    saveas(gcf,[plotsFolder,'geometric','oldVSnewAGQGdiffGPSpoints','.fig']) 
    saveas(gcf,[plotsFolder,'geometric','oldVSnewAGQGdiffGPSpoints','.png']) 
end