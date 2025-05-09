function computeParallelGravimetryLSC(GRID_PARA,COV_PARA,DEM_PARA,GRAV_PARA,OUTPUT_PARA,GRID_REF,Grav, ...
    GGM_Gravity_griddedInterpolant,ZDEM_griddedInterpolant,RTM_Correction_function, ...
    LWLBouguer_Slab_Function,density)
% computeGravimetryLSC runs the least squares collocation in blocks. 
% 
% Input:  GRID_PARA = grid/tiling parameters such as extent(MINLONG,MAXLONG,MINLAT,MAXLAT), buffer, STEP 
%         COV_PARA=   covariance function parameters
%         DEM_PARA=   DEM data such as filename,num_cols,num_rows,Regolith_filename?
%         GRAV_PARA=  gravity data such as filename,TypeB?,Grav_Faye_TypeB?
%         LEVELLING_PARA= levelling datacomparisons(Lev_eval,filename,Plot_Stat,Compare_To_Existing_Model,...)
%         OUTPUT_PARA= output parameters such as Grids_name,Tiles_dir_name
%         PLOT_GRIDS=True/False: gridded solution is plotted and output as well as the tiles
%         LongDEM= matrix of longitudes for DEM from Import_And_Format_Data_Sets_Reg
%         LatDEM=  matrix of latitudes for DEM from Import_And_Format_Data_Sets_Reg
%         GRID_REF= matrix of [*,3] from Import_And_Format_Data_Sets_Reg?
%         Grav = [longitude latitude orthometric_height gravity_anomaly uncertainty(std) flag(1:terrestrial 2:airborne 3:gradiometry)]
%         GGM_Gravity_griddedInterpolant = GGM gravity griddedInterpolant 
%         GGM_Zeta_griddedInterpolant = GGM height anomaly griddedInterpolant 
%         ZDEM_griddedInterpolant = DEM elevation griddedInterpolant
%         RTM_Correction_function=1*1 griddedInterpolant from RTM_CORRECTIONS?
%         LWLBouguer_Slab_Function=1*1 griddedInterpolant from RTM_CORRECTIONS?
%         Lev= ?from Import_And_Format_Data_Sets_Reg?
%         REFERENCE_GEOID_Zetai=1*1 griddedInterpolant from Import_And_Format_Data_Sets_Reg?
%         Coastline= 1*1 struct from Import_And_Format_Data_Sets_Reg
%
% Output:   function has no output.
% 
% Example: see Run_LSC_Drive_RTM
%
% Main functions
% - computeCovarianceFunctionParameters
% - precomputeCovarianceFunction
% - interpolateCovarianceFunction
% Other functions
% - 
% - 
%
% Written by Jack McCubbine
% Last updated by Neda Darbeheshti
% Geoscience Australia, 2024-03.
constants                                       % load constants
% Remove long wavelength of GGM
% Latitude (Grav(:,2)) has minus sign to be increased, height above the topography 
GGM_Gravity = GGM_Gravity_griddedInterpolant(Grav(:,1),-Grav(:,2),Grav(:,3)-ZDEM_griddedInterpolant(Grav(:,1),Grav(:,2)));

Grav(:,4) = Grav(:,4) - GGM_Gravity;

% plot LSC input
%plotCustomScatter(Grav(:,1),Grav(:,2),GGM_Gravity,GRID_PARA,2*GRID_PARA.buffer,'GGMgravity','mGal',OUTPUT_PARA.plotsFolder)

%plotCustomScatter(Grav(:,1),Grav(:,2),Grav(:,4),GRID_PARA,2*GRID_PARA.buffer,'LSCinputGravity','mGal',OUTPUT_PARA.plotsFolder)

% Run least squares collocation in blocks

lat_range = GRID_PARA.MAXLAT:-GRID_PARA.STEP:GRID_PARA.MINLAT;
long_range = GRID_PARA.MINLONG:GRID_PARA.STEP:GRID_PARA.MAXLONG;

n_lat = length(lat_range);
n_long = length(long_range);

% Total number of iterations
n_total = n_lat * n_long;

file_names = cell(n_total, 1); % Preallocate cell array for file names

parfor block_counter = 1:n_total
 %parfor block_counter = 1:4 % to test parfor
    BouguerConstant=0.0419; 
    EarthMajorAxis = 6378.137; 
    EarthMinorAxis = 6356.752; 
    AbsoluteGravityEquator_mgal = 9.7803267715*(10^5); 
    NormalGravityConstant = 0.001931851353; 
    EarthEccentricitySquared = 0.00669438002290; 
    % Convert linear index to 2D indices
    [lat_idx, long_idx] = ind2sub([n_lat, n_long], block_counter);
    
    % Get the actual latitude and longitude values
    LATsi = lat_range(lat_idx);
    LONGsi = long_range(long_idx);
    
    disp(['Working on block ',' centred at: ',num2str(LONGsi),', ',num2str(LATsi)])
   
    Nodes = struct('Latmin',LATsi,'Latmax',LATsi,'Longmin',LONGsi,'Longmax',LONGsi);

    % Step 1.1 Extract the grid points over the tile region
    disp('Extract the datasets over the tile region')
    CoordsMM=[Nodes.Longmin-GRID_PARA.buffer,Nodes.Latmin-GRID_PARA.buffer;...
              Nodes.Longmin-GRID_PARA.buffer,Nodes.Latmax+GRID_PARA.buffer;...
              Nodes.Longmax+GRID_PARA.buffer,Nodes.Latmax+GRID_PARA.buffer;...
              Nodes.Longmax+GRID_PARA.buffer,Nodes.Latmin-GRID_PARA.buffer;...
              Nodes.Longmin-GRID_PARA.buffer,Nodes.Latmin-GRID_PARA.buffer];
    
    INOUT=inpolygon(GRID_REF(:,1),GRID_REF(:,2),CoordsMM(:,1),CoordsMM(:,2)); % Whole tile zone mask - n.b. data on edges of tile are unreliable
    GRID_REF_dat=GRID_REF(INOUT==1,:);
    GRID_REF_dat(:,1)=round(GRID_REF_dat(:,1)*60)/60;
    GRID_REF_dat(:,2)=round(GRID_REF_dat(:,2)*60)/60;
    % Extract the data over the tile region
    CoordsMMinner=[Nodes.Longmin-GRID_PARA.buffer2,Nodes.Latmin-GRID_PARA.buffer2;...
                   Nodes.Longmin-GRID_PARA.buffer2,Nodes.Latmax+GRID_PARA.buffer2;...
                   Nodes.Longmax+GRID_PARA.buffer2,Nodes.Latmax+GRID_PARA.buffer2;...
                   Nodes.Longmax+GRID_PARA.buffer2,Nodes.Latmin-GRID_PARA.buffer2;...
                   Nodes.Longmin-GRID_PARA.buffer2,Nodes.Latmin-GRID_PARA.buffer2];

    INOUTinner=inpolygon(GRID_REF(:,1),GRID_REF(:,2),CoordsMMinner(:,1),CoordsMMinner(:,2)); % Inner zone where data are reliable mask
    
    disp('Extracting gravity residual signal in block')
    INOUT_grav=inpolygon(Grav(:,1),Grav(:,2),CoordsMM(:,1),CoordsMM(:,2));
    Gravdatout_m=Grav(INOUT_grav==1,:);% Extracts the full res data 
    % Assign to vairables to be used thorugh the rest of the code
    Gravdatout = struct('Long',Gravdatout_m(:,1),'Lat',Gravdatout_m(:,2),'ortho_H',Gravdatout_m(:,3),...
        'Grav_Anom',Gravdatout_m(:,4),'STD',Gravdatout_m(:,5));

    % Step 1.2 Covariance analysis
    COV_PARA_RTM=COV_PARA;
    COV_PARA_Faye=COV_PARA;
    if COV_PARA.COV_COMPUTED_Tilewise==true
   
    Gravdatout_m_cv = Gravdatout_m;
   
    disp('Running localised covariance analysis') 

    if COV_PARA.Airbornedataonly && isempty(Gravdatout_m(Gravdatout_m(:,6)==2,:))==false

    [COV_PARA_RTM]= computeCovarianceFunctionParameters(OUTPUT_PARA,COV_PARA,Gravdatout_m_cv(Gravdatout_m(:,6)==2,:),'RTM',block_counter);
    
    Gravdatout_m_cv_faye=Gravdatout_m_cv;
    
    Gravdatout_m_cv_faye(:,4)=Gravdatout_m_cv(Gravdatout_m(:,6)==2,4)-LWLBouguer_Slab_Function(Gravdatout_m_cv(Gravdatout_m(:,6)==2,1),Gravdatout_m_cv(Gravdatout_m(:,6)==2,2))+...
        ZDEM_griddedInterpolant(Gravdatout_m_cv(Gravdatout_m(:,6)==2,1),Gravdatout_m_cv(Gravdatout_m(:,6)==2,2))*BouguerConstant*density;
    
    [COV_PARA_Faye]= computeCovarianceFunctionParameters(OUTPUT_PARA,COV_PARA,Gravdatout_m_cv_faye,'Faye',block_counter);
    
    else
    
        [COV_PARA_RTM] =  computeCovarianceFunctionParameters(OUTPUT_PARA,COV_PARA,Gravdatout_m_cv,'RTM',block_counter);    
       
        Gravdatout_m_cv_faye = Gravdatout_m_cv;
    
        Gravdatout_m_cv_faye(:,4)=Gravdatout_m_cv(:,4)-LWLBouguer_Slab_Function(Gravdatout_m_cv(:,1),Gravdatout_m_cv(:,2))+...
                ZDEM_griddedInterpolant(Gravdatout_m_cv(:,1),Gravdatout_m_cv(:,2))*BouguerConstant*density;
    
        [COV_PARA_Faye] =  computeCovarianceFunctionParameters(OUTPUT_PARA,COV_PARA,Gravdatout_m_cv_faye,'Faye',block_counter);    
    end
    end

    % radius of the Bjerhammar sphere
    phi=deg2rad(mean(GRID_REF_dat(:,2)));
    RadiusBjerhammar= EarthMajorAxis*EarthMinorAxis/sqrt((EarthMajorAxis*sin(phi)).^2+(EarthMinorAxis*cos(phi)).^2)*10^3;% Pajama sphere radius.
    
    disp('Pre-computing RTM gravimetry covariance function') 
    CCov_tt_int_fun_RTM=precomputeCovarianceFunction('cov_tt',RadiusBjerhammar,COV_PARA_RTM.width,COV_PARA_RTM.res,COV_PARA_RTM.A,COV_PARA_RTM.B,COV_PARA_RTM.N,COV_PARA_RTM.M);
    CCov_gg_int_fun_RTM=precomputeCovarianceFunction('cov_gg',RadiusBjerhammar,COV_PARA_RTM.width,COV_PARA_RTM.res,COV_PARA_RTM.A,COV_PARA_RTM.B,COV_PARA_RTM.N,COV_PARA_RTM.M);
    CCov_gt_int_fun_RTM=precomputeCovarianceFunction('cov_gt',RadiusBjerhammar,COV_PARA_RTM.width,COV_PARA_RTM.res,COV_PARA_RTM.A,COV_PARA_RTM.B,COV_PARA_RTM.N,COV_PARA_RTM.M);
      
    disp('Pre-computing Faye gravimetry covariance function') 
    CCov_gg_int_fun_Faye=precomputeCovarianceFunction('cov_gg',RadiusBjerhammar,COV_PARA_Faye.width,COV_PARA_Faye.res,COV_PARA_Faye.A,COV_PARA_Faye.B,COV_PARA_Faye.N,COV_PARA_Faye.M);
    CCov_gt_int_fun_Faye=precomputeCovarianceFunction('cov_gt',RadiusBjerhammar,COV_PARA_Faye.width,COV_PARA_Faye.res,COV_PARA_Faye.A,COV_PARA_Faye.B,COV_PARA_Faye.N,COV_PARA_Faye.M);
        
    % Compute normal gravity at the computation points on the tile

    NormalGravity=AbsoluteGravityEquator_mgal*(1+NormalGravityConstant*(sin(deg2rad(GRID_REF_dat(:,2))).^2) ...
        )./sqrt(1-EarthEccentricitySquared*(sin(deg2rad(GRID_REF_dat(:,2))).^2));
   
    disp(['The size of the LSC matrix is going to be: ',num2str(length(Gravdatout.Long)),...
        ' x ',num2str(length(Gravdatout.Long))])
   
    if ~isempty(Gravdatout.Long)
        disp('... Interpolating gravimetry covariance matrices')
      
        % Auto-covariance of gravity anomaly at gravity points
        ACOVggRTM_grav = interpolateCovarianceFunction(...
            Gravdatout.Long, Gravdatout.Lat, RadiusBjerhammar + Gravdatout.ortho_H, ...
            Gravdatout.Long, Gravdatout.Lat, RadiusBjerhammar + Gravdatout.ortho_H, CCov_gg_int_fun_RTM);
    
        % Auto-covariance of gravity anomaly at DEM points
        ACOVggRTM_dem = interpolateCovarianceFunction(...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), ...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), CCov_gg_int_fun_RTM);
    
        % Cross-covariance of gravity anomaly at gravity and DEM points
        CCOVggRTM_dem_grav = interpolateCovarianceFunction(...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), ...
            Gravdatout.Long, Gravdatout.Lat, RadiusBjerhammar + Gravdatout.ortho_H, CCov_gg_int_fun_RTM);
    
        % Cross-covariance of gravity anomaly and potential at gravity and DEM points
        CCOVgtRTM_dem_grav = interpolateCovarianceFunction(...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), ...
            Gravdatout.Long, Gravdatout.Lat, RadiusBjerhammar + Gravdatout.ortho_H, CCov_gt_int_fun_RTM);
        
        % Auto-covariance of potential at DEM points
        ACOVttRTM_dem = interpolateCovarianceFunction(...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), ...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), CCov_tt_int_fun_RTM);
        
        % Auto-covariance of Faye gravity anomaly at DEM points
        ACOVggFay_dem = interpolateCovarianceFunction(...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), ...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), CCov_gg_int_fun_Faye);
    
        % Cross-covariance of Faye gravity anomaly and potential at DEM points
        CCOVgtFay_dem = interpolateCovarianceFunction(...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), ...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), CCov_gt_int_fun_Faye);

        disp('... Solving LSC matrix equations for Gravimetry') 

        [LSCResidualGeoid,errorMatrixGeoid,residualFreeAirGravityAnomaly,residualBouguerGravityAnomaly,...
               errorMatrixGravity]=solveGravityLSCmatrix(GRAV_PARA,Gravdatout,GRID_REF_dat,density, ...
               ACOVggRTM_grav,ACOVggRTM_dem,CCOVggRTM_dem_grav,CCOVgtRTM_dem_grav,ACOVttRTM_dem,ACOVggFay_dem, ...
               CCOVgtFay_dem,RTM_Correction_function,LWLBouguer_Slab_Function,ZDEM_griddedInterpolant,NormalGravity);

        disp('Done')    

       % Create weights for grid blending 
    
        inputWeightVector=INOUT.*INOUTinner;
    
        filterWeights = createGridWeights(inputWeightVector, DEM_PARA, GRID_PARA);
    
        % Save the data to a tile 
    
        weights=filterWeights;
    
        res_geoid=INOUT*0;
        res_geoid(INOUT==1)=LSCResidualGeoid;
    
        pot_error=INOUT*0;
        pot_error(INOUT==1)=abs(errorMatrixGeoid);
        
        res_grav=INOUT*0;
        res_grav(INOUT==1)=residualFreeAirGravityAnomaly;
    
        res_grav_Bouguer=INOUT*0;
        res_grav_Bouguer(INOUT==1)=residualBouguerGravityAnomaly;
    
        grav_error=INOUT*0;
        grav_error(INOUT==1)=abs(errorMatrixGravity);
    
        Dataset_save = struct('weights',weights,'res_geoid',res_geoid,'pot_error',pot_error,...
        'res_grav',res_grav,'res_grav_Bouguer',res_grav_Bouguer,'grav_error',grav_error,...
        'COV_PARA_RTM', COV_PARA_RTM,'COV_PARA_Faye',COV_PARA_Faye);

        % Construct the file path
        file_names{block_counter} = fullfile(OUTPUT_PARA.Tiles_dir_name, ['Tile', num2str(LONGsi), '_', num2str(abs(LATsi)), '.mat']);
     
        parsave(file_names{block_counter},Dataset_save);
    
    else
    disp(('No data in block'))
    end
end
%
end