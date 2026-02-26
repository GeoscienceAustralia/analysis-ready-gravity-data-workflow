function customizeMap(quantityName,colorbarUnit,Coastline,axisLimits)
% Function to customize each plot
   
    for jj = 1:9
        plot(Coastline.long{jj}, Coastline.lat{jj}, 'k')
    end

    colorbar;
    colormap(gmt_haxby(256));
    title(colorbar,colorbarUnit,'FontSize',10);
    axis([axisLimits.lonMinLimit axisLimits.lonMaxLimit axisLimits.latMinLimit axisLimits.latMaxLimit])
    ax = gca;
    ax.PlotBoxAspectRatio =[1 axisLimits.latMeanCosine 1];
    title(quantityName)
    xlabel('Longitude')
    ylabel('Latitude')
   
end