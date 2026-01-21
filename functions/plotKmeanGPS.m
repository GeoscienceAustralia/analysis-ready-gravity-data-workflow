function plotKmeanGPS(GPSlevelling3D,geomGravDiff,geomRefAGQGDiff,Coastline,GRID_PARA,plotsFolder)

    % common variables for plotting
    axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
    axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
    axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
    axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
    axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;

    % plot GPSlevelling 
    figure('Name','MosaicTiles','NumberTitle','off'); 
    clf
    hold on
    scatter(GPSlevelling3D(:,1),GPSlevelling3D(:,2),15,GPSlevelling3D(:,3),'filled')
    customizeMap('GPSlevelling','m',Coastline,axisLimits)
    saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevelling','.png']) 

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
    customizeMap('Geometric and LSC AGQG Difference','m',Coastline,axisLimits)
    caxis([-0.3 0.3])
    saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevellingLSC','.png']) 
    
    % Plot GPS levelling vs LSC with clustering
%     figure('Name','MosaicTiles','NumberTitle','off'); 
%     clf
%     hold on
%     scatter(GPSlevelling3D(:,1), GPSlevelling3D(:,2), 15, idx, 'filled');
%     customizeMap('Clustered GPS levelling points',' ',Coastline,axisLimits)
%     colorbar off;
%     saveas(gcf,[plotsFolder,'MosaicTiles','clusteredGPSlevellingLSC','.png']) 
       
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
    customizeMap('Geometric and 2022 AGQG Difference','m',Coastline,axisLimits)
    caxis([-0.3 0.3])
    saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevellingAGQG','.png']) 
    
    % Plot GPS levelling vs reference AGQG with clustering
%     figure('Name','MosaicTiles','NumberTitle','off'); 
%     clf
%     hold on
%     scatter(GPSlevelling3D(:,1), GPSlevelling3D(:,2), 15, idx, 'filled');
%     customizeMap('Clustered GPS levelling points',' ',Coastline,axisLimits)
%     colorbar off;
%     saveas(gcf,[plotsFolder,'MosaicTiles','clusteredGPSlevellingAGQG','.png']) 


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

% Remove a tiled plane so the signal is zero mean for the LSC

% Construct the matrix for linear trend removal
trendMatrix = [GPSlevelling3D(:,1) - mean(GPSlevelling3D(:,1)), GPSlevelling3D(:,2) - mean(GPSlevelling3D(:,2)), ones(size(GPSlevelling3D(:,2)))];

% Calculate the coefficients of the best-fit plane
trendCoefficients = trendMatrix \ geomGravDiff;
trendCoefficients2022 = trendMatrix \ geomRefAGQGDiff;

% Remove the planar trend to obtain zero-mean data
geomGravDiffDetrended = geomGravDiff - trendMatrix * trendCoefficients;
geomGravDiffDetrended2022 = geomRefAGQGDiff - trendMatrix * trendCoefficients2022;

% plot differences between geometric and gravimetric geoid at GPS leveling points 

figure('Name','MosaicTiles','NumberTitle','off'); 
clf
hold on
scatter(GPSlevelling3D(:,1),GPSlevelling3D(:,2),15,geomGravDiffDetrended,'filled')
customizeMap('Geometric and LSC AGQG Difference Detrended','m',Coastline,axisLimits)
caxis([-0.3 0.3])
saveas(gcf,[plotsFolder,'MosaicTiles','geomGravDiffDetrended','.png']) 

figure('Name','MosaicTiles','NumberTitle','off'); 
clf
hold on
scatter(GPSlevelling3D(:,1),GPSlevelling3D(:,2),15,geomGravDiffDetrended2022,'filled')
customizeMap('Geometric and 2022 AGQG Difference Detrended','m',Coastline,axisLimits)
caxis([-0.3 0.3])
saveas(gcf,[plotsFolder,'MosaicTiles','geomGravDiffDetrended2022','.png']) 

figure('Name','MosaicTiles','NumberTitle','off'); 
clf
hold on
scatter(GPSlevelling3D(:,1),GPSlevelling3D(:,2),15,geomGravDiff,'filled')
customizeMap('Geometric and LSC AGQG Difference','m',Coastline,axisLimits)
saveas(gcf,[plotsFolder,'MosaicTiles','geomGravDiff','.png']) 

figure('Name','MosaicTiles','NumberTitle','off'); 
clf
hold on
scatter(GPSlevelling3D(:,1),GPSlevelling3D(:,2),15,geomRefAGQGDiff,'filled')
customizeMap('Geometric and 2022 AGQG Difference','m',Coastline,axisLimits)
saveas(gcf,[plotsFolder,'MosaicTiles','geomGravDiff','.png'])

fprintf('%f length  Geometric and LSC AGQG Difference Detrended\n',length (geomGravDiffDetrended));
fprintf('%f min     Geometric and LSC AGQG Difference Detrended\n',min    (geomGravDiffDetrended));
fprintf('%f max     Geometric and LSC AGQG Difference Detrended\n',max    (geomGravDiffDetrended));
fprintf('%f mean    Geometric and LSC AGQG Difference Detrended\n',mean   (geomGravDiffDetrended));
fprintf('%f median  Geometric and LSC AGQG Difference Detrended\n',median (geomGravDiffDetrended));
fprintf('%f std     Geometric and LSC AGQG Difference Detrended\n',std    (geomGravDiffDetrended));


fprintf('%f length  Geometric and 2022 AGQG Difference Detrended\n',length (geomGravDiffDetrended2022));
fprintf('%f min     Geometric and 2022 AGQG Difference Detrended\n',min    (geomGravDiffDetrended2022));
fprintf('%f max     Geometric and 2022 AGQG Difference Detrended\n',max    (geomGravDiffDetrended2022));
fprintf('%f mean    Geometric and 2022 AGQG Difference Detrended\n',mean   (geomGravDiffDetrended2022));
fprintf('%f median  Geometric and 2022 AGQG Difference Detrended\n',median (geomGravDiffDetrended2022));
fprintf('%f std     Geometric and 2022 AGQG Difference Detrended\n',std    (geomGravDiffDetrended2022));








