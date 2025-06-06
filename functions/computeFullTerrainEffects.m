function [fullTopographyCorrected_gravityPoint,longwaveTopography_griddedInterpolant,fullTopography_griddedInterpolant,fullTopographyCorrected_gravityGradient]=computeFullTerrainEffects(GRID_PARA, ...
    Topo_PARA,Grav,Grav_grad,GGM_Gravity_griddedInterpolant,DEM_data,ZDEM_griddedInterpolant,LongDEM_matrix,LatDEM_matrix,Coastline,plotsFolder)
    % computeTerrainEffect calculates the residual terrain model. It first calculates
    % the full Bouguer correction, then filters the DEM and minimises the
    % residual gravity signal w.r.t to the GGM.

    % Input:  GRID_PARA = struct with buffer,buffer2,STEP,MINLONG,MAXLONG,MINLAT,MAXLAT 
    %         Topo_PARA = struct with Corr, TopoPlot,Density,Depth,Rad,RTM[,,] 
    %         GRAV = [latitude longitude height gravity uncertainty flag]
    %         GGM_Gravity_griddedInterpolant = GGM gravity griddedInterpolant 
    %         DEM_data = [longitudeDEM latitudeDEM heightDEM] 
    %         ZDEM_griddedInterpolant = DEM elevation griddedInterpolant
    %         LongDEM_matrix = matrix of DEM longitude [DEM_PARA.num_rows*DEM_PARA.num_cols] 
    %         LatDEM_matrix = matrix of DEM latitude [3181*4861] 
 
    % Output:   fullTopographyCorrected_gravityPoint = [latitude longitude height gravity_rtm uncertainty flag]
    %           longwaveTopography_griddedInterpolant =  1*1 griddedInterpolant
    %           fullTopography_griddedInterpolant = 1*1 griddedInterpolant

    % Example: see Run_LSC_Drive_RTM

    % Main functions
    % - computeTerrainCorrection: computes the terrain correction.
    % - filterDEM: 
    % Other functions
    % - inpolygon:  
    % - std: 
    % - unique: 

    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11.

    constants                                       % load constants
    
    LatDEM_topo  = reshape(DEM_data(:,2),length(unique(DEM_data(:,1))),length(unique(DEM_data(:,2))))';
    LongDEM_topo = reshape(DEM_data(:,1),length(unique(DEM_data(:,1))),length(unique(DEM_data(:,2))))';
    ZDEM_topo    = reshape(DEM_data(:,3),length(unique(DEM_data(:,1))),length(unique(DEM_data(:,2))))';

    disp('Running topographic corrections for gravity points')
    TC_gravity_point = computeTerrainCorrection(ZDEM_topo,LatDEM_topo,LongDEM_topo,Grav(:,2),Grav(:,1),Grav(:,3),     Topo_PARA.Rad,Topo_PARA.Density,'g');
    %save([plotsFolder, 'terrainEffects1','.mat'], 'TC_gravity_point'); 
    disp('Running topographic corrections for DEM points')
    TC_DEM_point =     computeTerrainCorrection(ZDEM_topo,LatDEM_topo,LongDEM_topo,LatDEM_topo,LongDEM_topo,ZDEM_topo,Topo_PARA.Rad,Topo_PARA.Density,'g');
    %save([plotsFolder, 'terrainEffects2','.mat'], 'TC_DEM_point');
    if Topo_PARA.TopoPlot
        
%         plotProfiles(Grav(13676:13796,2),Grav(13676:13796,3),TC_gravity_point(13676:13796,1),Grav(13676,1),'Latitude','Height [m]','Gravity correction [mGal]','PrismGravityEffectHeight',plotsFolder)
%         plotProfiles(LatDEM_topo(:,90),ZDEM_topo(:,90),TC_DEM_point(:,90),LongDEM_topo(1,90),'Latitude','DEM[m]','GravityCorrection[mGal]','PrismGravityEffectDEM',plotsFolder)
        plotCustomScatter(LongDEM_topo, LatDEM_topo, TC_DEM_point, GRID_PARA,'DEMPrismGravityEffect','mGal',Coastline,[],plotsFolder);
        plotCustomScatter(Grav(:,1),Grav(:,2),TC_gravity_point, GRID_PARA,'gravityTopoPrismGravityEffect','mGal',Coastline,[],plotsFolder);
      
    end
    
    disp('Optimising DEM filter for spherical harmonics degree ')

    % Estimate the maximum number of iterations based on the range of spherical harmonic degrees
    max_iterations = numel(Topo_PARA.RTM(1):Topo_PARA.RTM(2):Topo_PARA.RTM(3));
    
    % Preallocate Minimise_WL with two zero columns (N and std(gravity_residual))
    Minimise_WL = zeros(max_iterations, 2); 
    
    % DEM latitude and longitude extension
    minLatDEM = min(min(LatDEM_matrix));
    maxLatDEM = max(max(LatDEM_matrix));
    minLonDEM = min(min(LongDEM_matrix));
    maxLonDEM = max(max(LongDEM_matrix));
    
    latExtent = [minLatDEM, maxLatDEM];
    lonExtent = [minLonDEM, maxLonDEM];
    
    % Loop over different spherical harmonic degrees
    for iteration_counter = 1:max_iterations
        N = Topo_PARA.RTM(1) + (iteration_counter - 1) * Topo_PARA.RTM(2);
        
        % Filter DEM and calculate LWLBouguer_Slab_Function
        [filtered_dem] = filterDEM(ZDEM_griddedInterpolant(LongDEM_matrix, LatDEM_matrix), N, latExtent, lonExtent);
        longwaveTopography_griddedInterpolant = griddedInterpolant(LongDEM_matrix', LatDEM_matrix(end:-1:1, :)', BouguerConstant * Topo_PARA.Density * filtered_dem(end:-1:1, :)');
        
        % Calculate gravity residuals
        gravity_residual = (Grav(:, 4) - TC_gravity_point) - (GGM_Gravity_griddedInterpolant(Grav(:, 1), -Grav(:, 2), Grav(:, 3)) ...
            - longwaveTopography_griddedInterpolant(Grav(:, 1), Grav(:, 2)));
        
        % Append the harmonic degree (N) and the standard deviation of gravity residuals to Minimise_WL
        Minimise_WL(iteration_counter, :) = [N, std(gravity_residual)];
    end
    % outputs
    % A grid function of Bouguer slab corrections
    [~,i]=min(Minimise_WL(:,2));
    N_min=Minimise_WL(i,1);

    disp(N_min)
    
    [filtered_dem] = filterDEM(ZDEM_griddedInterpolant(LongDEM_matrix,LatDEM_matrix), N_min, latExtent, lonExtent);
    longwaveTopography_griddedInterpolant=griddedInterpolant(LongDEM_matrix',LatDEM_matrix(end:-1:1,:)',BouguerConstant*Topo_PARA.Density*filtered_dem(end:-1:1,:)');
    
    % Subtract terrain and Bouguer slab corrections from gravity
    fullTopographyCorrected_gravityPoint=Grav;
    fullTopographyCorrected_gravityPoint(:,4)=fullTopographyCorrected_gravityPoint(:,4)-TC_gravity_point+longwaveTopography_griddedInterpolant(Grav(:,1),Grav(:,2));
     
    % A grid function of sum of terrain and Bouguer slab corrections at DEM
    fullTopographyCorrection_dem=-TC_DEM_point+longwaveTopography_griddedInterpolant(LongDEM_topo,LatDEM_topo);
    fullTopography_griddedInterpolant=griddedInterpolant(LongDEM_topo',LatDEM_topo(end:-1:1,:)',fullTopographyCorrection_dem(end:-1:1,:)');
    
  
    fullTopographyCorrected_gravityGradient=Grav_grad;

    if ~isempty(Grav_grad)
        
       disp('Running topographic corrections for gravity gradient')

       % Compute gravity gradient corrections at gradiometry points
       TCgradient_gradient_point = computeTerrainCorrection(ZDEM_topo,LatDEM_topo,LongDEM_topo,Grav_grad(:,2),Grav_grad(:,1),Grav_grad(:,3),     Topo_PARA.Rad,Topo_PARA.Density,'gg');
      
       % plot gravity gradient prism topography correction
       % plotProfiles(Grav(13676:13796,2),Grav(13676:13796,3),-TCgradient_gravity_point(13676:13796,1),Grav(13676,1),'Latitude','Height [m]','Gravity gradient correction [mGal/m]','PrismGravityGradientEffectHeight',plotsFolder)
       plotCustomScatter(Grav_grad(:,1),Grav_grad(:,2),TCgradient_gradient_point, GRID_PARA,'flightAltitudePrismGravityGradientEffect','mGal/m',Coastline,[],plotsFolder);
    
       % Subtract gravity gradient prism topography correction from gravity gradient
       fullTopographyCorrected_gravityGradient(:,4)=fullTopographyCorrected_gravityGradient(:,4)-TCgradient_gradient_point;
    end

end