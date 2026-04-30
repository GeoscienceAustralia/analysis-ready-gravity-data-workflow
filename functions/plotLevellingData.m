function plotLevellingData(LonLat,Lev,quantityName,colorbarUnit,plotsFolder)

    figure('Name','GPSlevelling','NumberTitle','off'); 
    clf
    hold on
    scatter(LonLat(:,1),LonLat(:,2),5,Lev,'filled')
    colormap(jet);
    %caxis([0 1])
    title(colorbar,colorbarUnit,'FontSize',10);
    colorbar;
    title(quantityName)
    xlabel('Longitude')
    ylabel('Latitude')
    saveas(gcf,[plotsFolder,'GPSlevelling',titlePlot,'.fig']) 
    saveas(gcf,[plotsFolder,'GPSlevelling',titlePlot,'.png'])

 end