function gravity6D = weightInputGravity(gravity6D,Coastline,GRID_PARA,OUTPUT_PARA)

    % plot gravity data with different source
    
    gravFlagTerrestrial  = gravity6D(:,6)==1;
    gravFlagAltimetry  = gravity6D(:,6)==2;
    gravFlagAirborne  = gravity6D(:,6)==3;

    if OUTPUT_PARA.PLOT_GRIDS

        plotCustomScatter(gravity6D(gravFlagTerrestrial,1),gravity6D(gravFlagTerrestrial,2),gravity6D(gravFlagTerrestrial,5),GRID_PARA,'Terrestrialuncertainty','mGal',Coastline,[],OUTPUT_PARA.plotsFolder)
        plotCustomScatter(gravity6D(gravFlagAltimetry,1),gravity6D(gravFlagAltimetry,2),gravity6D(gravFlagAltimetry,5),GRID_PARA,'Altimetryuncertainty','mGal',Coastline,[],OUTPUT_PARA.plotsFolder)
        plotCustomScatter(gravity6D(gravFlagAirborne,1),gravity6D(gravFlagAirborne,2),gravity6D(gravFlagAirborne,5),GRID_PARA,'Airborneuncertainty','mGal',Coastline,[],OUTPUT_PARA.plotsFolder)

        plotCustomScatter(gravity6D(gravFlagTerrestrial,1),gravity6D(gravFlagTerrestrial,2),gravity6D(gravFlagTerrestrial,3),GRID_PARA,'TerrestrialHeight','m',Coastline,[],OUTPUT_PARA.plotsFolder)
        plotCustomScatter(gravity6D(gravFlagAltimetry,1),gravity6D(gravFlagAltimetry,2),gravity6D(gravFlagAltimetry,3),GRID_PARA,'AltimetryHeight','m',Coastline,[],OUTPUT_PARA.plotsFolder)
        plotCustomScatter(gravity6D(gravFlagAirborne,1),gravity6D(gravFlagAirborne,2),gravity6D(gravFlagAirborne,3),GRID_PARA,'AirborneHeight','m',Coastline,[],OUTPUT_PARA.plotsFolder)

    end
    % mark airborne in altimetry
    % Extract data inside computational extents
    % airborne.MINLONG=153.3;
    % airborne.MAXLONG=153.8;
    % airborne.MINLAT=-28.1;
    % airborne.MAXLAT=-30;
    % disp('Extracting airborne gravity subset') 
    % CoordsMM_Grav=[airborne.MINLONG,airborne.MINLAT;...
    %                airborne.MINLONG,airborne.MAXLAT;...
    %                airborne.MAXLONG,airborne.MAXLAT;...
    %                airborne.MAXLONG,airborne.MINLAT;...
    %                airborne.MINLONG,airborne.MINLAT];
    %         
    % Gravin=inpolygon(Gravo(:,1),Gravo(:,2),CoordsMM_Grav(:,1),CoordsMM_Grav(:,2));
    % % case 1: Eliminate rows where the 6th column equals 2 and Gravin is true
    % % Gravo = Gravo(~(Gravo(:,6) == 2 & Gravin), :);
    % % case 2: Add more uncertainty to the fifth (uncertainty) column where altimetry data 
    %   Gravo( Gravo(:,6) == 2 & Gravin, 5) = Gravo( Gravo(:,6) == 2 & Gravin , 5) + 100;
    % 
    % plotCustomScatter(Gravo(:,1),Gravo(:,2),Gravo(:,5),GRID_PARA,'GravityUncertainty','mGal',Coastline,[],OUTPUT_PARA.plotsFolder)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% eliminate data based on data flags
    
    gravity6D(gravFlagAltimetry==1 | gravFlagTerrestrial==1,:)=[]; % for two data types
    
     %gravity6D(gravFlagAirborne==1,:)=[]; %just one data type

    if OUTPUT_PARA.PLOT_GRIDS
        plotCustomScatter(gravity6D(:,1),gravity6D(:,2),gravity6D(:,5),GRID_PARA,'GravityUncertainty','mGal',Coastline,[],OUTPUT_PARA.plotsFolder)
        plotCustomScatter(gravity6D(:,1),gravity6D(:,2),gravity6D(:,4),GRID_PARA,'Gravity','mGal',Coastline,[],OUTPUT_PARA.plotsFolder)
    end
end