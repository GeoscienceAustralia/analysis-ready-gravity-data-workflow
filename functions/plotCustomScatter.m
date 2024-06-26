function plotCustomScatter(longVector, latvector, dataVector, GRID_PARA, GridRadius,quantityName,dataunit,plotsFolder)
    % Input:  longVector   = vector of longitudes
    %         latvector    = vector of latitudes
    %         dataVector   = vector of data
    %         GRID_PARA    = grid parameters such as extent (MINLONG,MAXLONG,MINLAT,MAXLAT) and buffer
    %         GridRadius   = number in degrees
    %         quantityName = character
    %         dataunit     = character
    %         plotFolder   = character
    % 
    % Output: saved plot in png format

    % Example: 
    % plotCustomScatter(DEM_data(:,1),DEM_data(:,2),DEM_data(:,3),GRID_PARA,Topo_PARA.Rad,'DEM','m',plotsFolder)
    %
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2024-02.
  
    % plot scatter data
                
    figure('Name','plotScatter','NumberTitle','off');
    clf
    hold on
    scatter(longVector(:),latvector(:),1,dataVector(:))
    colormap(jet)
    colorbar
    title(colorbar,dataunit,'FontSize',10);
    axis([GRID_PARA.MINLONG-2*GridRadius GRID_PARA.MAXLONG+2*GridRadius GRID_PARA.MINLAT-2*GridRadius GRID_PARA.MAXLAT+2*GridRadius])
    ax=gca;
    ax.PlotBoxAspectRatio =[1 abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT])))) 1];
    title(quantityName)
    xlabel('Longitude')
    ylabel('Latitude')
    saveas(gcf,[plotsFolder,'scatter',quantityName,'.png'])
end