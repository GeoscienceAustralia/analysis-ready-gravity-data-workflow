% Plot GPS levelling vs LSC with clustering
figure('Name','MosaicTiles','NumberTitle','off'); 
clf
hold on

% Clean and center the data
validVals = geomGravDiff(~isnan(geomGravDiff));
valDiff = geomGravDiff - mean(validVals);

% Apply K-means clustering
k = 3;
[idx, C] = kmeans(valDiff, k);

% Plot clustered data
scatter(Lev(:,1), Lev(:,2), 15, idx, 'filled');
colormap(lines(k)); % Use distinguishable colors
caxis([1 k]); % Match color axis to number of clusters

% Remove colorbar (scale bar)
colorbar off;

% Add legend
legend(arrayfun(@(x) sprintf('Cluster %d', x), 1:k, 'UniformOutput', false));

title('Clustered GPS Levelling vs LSC');
xlabel('Longitude');
ylabel('Latitude');
hold off;



% Initialize statistics containers
for i = 1:k
    clusterVals = valDiff(idx == i);
    fprintf('Cluster %d Statistics:\n', i);
    fprintf('  Count: %d\n', numel(clusterVals));
    fprintf('  Mean: %.4f\n', mean(clusterVals));
    fprintf('  Std Dev: %.4f\n', std(clusterVals));
    fprintf('  Min: %.4f\n', min(clusterVals));
    fprintf('  Max: %.4f\n\n', max(clusterVals));
end








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
%customizeMap('Geometric and AGQG Difference','m',Coastline,axisLimits)
caxis([-0.1 0.1])
%saveas(gcf,[plotsFolder,'MosaicTiles','GPSlevellingAGQG','.fig']) 

validVals = AGQG_Vals_Lev(~isnan(AGQG_Vals_Lev));
valDiff = AGQG_Vals_Lev - mean(validVals);

fprintf('%f length  GPSlevellingAGQG\n',length (valDiff));
fprintf('%f min     GPSlevellingAGQG\n',min    (valDiff));
fprintf('%f max     GPSlevellingAGQG\n',max    (valDiff));
fprintf('%f mean    GPSlevellingAGQG\n',mean   (valDiff));
fprintf('%f median  GPSlevellingAGQG\n',median (valDiff));
fprintf('%f std     GPSlevellingAGQG\n',std    (valDiff));
