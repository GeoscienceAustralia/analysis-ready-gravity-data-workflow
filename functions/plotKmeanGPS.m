function plotKmeanGPS(GPSlevelling3D,geomGravDiff,geomRefAGQGDiff,Coastline,GRID_PARA,plotsFolder)

    % common variables for plotting
    axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
    axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
    axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
    axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
    axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;

    % Clean and center the data
    validVals = geomGravDiff(~isnan(geomGravDiff));
    valDiff = geomGravDiff - mean(validVals);

    % Apply K-means clustering
    k = 1;
    [idx, ~] = kmeans(valDiff, k);

    % plot GPSlevelling vs LSC
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    scatter(GPSlevelling3D(:,1),GPSlevelling3D(:,2),15,valDiff,'filled')
    customizeMap('Geometric and LSC Gravimetric Geoid Difference','m',Coastline,axisLimits)
    saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevellingLSC','.png']) 
    
    % Plot GPS levelling vs LSC with clustering
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    scatter(GPSlevelling3D(:,1), GPSlevelling3D(:,2), 15, idx, 'filled');
    customizeMap('Clustered GPS levelling points',' ',Coastline,axisLimits)
    colorbar off;
    saveas(gcf,[plotsFolder,'MosaicTiles','clusteredGPSlevellingLSC','.png']) 
       
    % Initialize statistics containers
    for i = 1:k
        clusterVals = valDiff(idx == i);
        fprintf('Cluster Geometric and LSC Gravimetric Geoid Difference %d Statistics:\n', i);
        fprintf('  Count: %d\n', numel(clusterVals));
        fprintf('  Mean: %.4f\n', mean(clusterVals));
        fprintf('  Std Dev: %.4f\n', std(clusterVals));
        fprintf('  Min: %.4f\n', min(clusterVals));
        fprintf('  Max: %.4f\n\n', max(clusterVals));
    end


    % Clean and center the data
    validVals = geomRefAGQGDiff(~isnan(geomRefAGQGDiff));
    valDiff = geomRefAGQGDiff - mean(validVals);

    % plot GPSlevelling vs reference AGQG
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    scatter(GPSlevelling3D(:,1),GPSlevelling3D(:,2),15,valDiff,'filled')
    customizeMap('Geometric and AGQG Difference','m',Coastline,axisLimits)
   saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevellingAGQG','.png']) 
    
    % Plot GPS levelling vs reference AGQG with clustering
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    scatter(GPSlevelling3D(:,1), GPSlevelling3D(:,2), 15, idx, 'filled');
    customizeMap('Clustered GPS levelling points',' ',Coastline,axisLimits)
    colorbar off;
    saveas(gcf,[plotsFolder,'MosaicTiles','clusteredGPSlevellingAGQG','.png']) 


    % Initialize statistics containers
    for i = 1:k
        clusterVals = valDiff(idx == i);
        fprintf('Cluster Geometric and AGQG Difference %d Statistics:\n', i);
        fprintf('  Count: %d\n', numel(clusterVals));
        fprintf('  Mean: %.4f\n', mean(clusterVals));
        fprintf('  Std Dev: %.4f\n', std(clusterVals));
        fprintf('  Min: %.4f\n', min(clusterVals));
        fprintf('  Max: %.4f\n\n', max(clusterVals));
    end


