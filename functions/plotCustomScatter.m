function plotCustomScatter(longVector, latvector, dataVector, GRID_PARA, quantityName, dataunit, Coastline, plotsFolder)
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

    % common variables for plotting
    axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
    axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
    axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
    axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
    axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;
  
    % plot scatter data
                
    figure('Name','plotScatter','NumberTitle','off');
    clf
    hold on
    scatter(longVector(:),latvector(:),.9,dataVector(:))
    customizeMap(quantityName,dataunit,Coastline,axisLimits)
    caxis([-100 150])
    saveas(gcf,[plotsFolder,'scatter',quantityName,'.png'])
end