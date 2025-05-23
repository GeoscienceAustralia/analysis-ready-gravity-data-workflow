function computeGradiometryLSC(GRID_PARA,COV_PARA,DEM_PARA,GRAV_PARA,GRAV_GRAD_PARA,OUTPUT_PARA,GRID_REF,Grav,Grav_grad, ...
    GGM_Gravity_griddedInterpolant,ZDEM_griddedInterpolant,RTM_Correction_function, ...
    LWLBouguer_Slab_Function,density)
% computeGravimetryGradiometryLSC runs the least squares collocation in blocks. 
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
%         Grav_grad= [longitude latitude height gravity_gradient flag (5:gradiometry)]
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

if GRAV_GRAD_PARA.avail
    GGM_GravityGradient=(GGM_Gravity_griddedInterpolant(Grav_grad(:,1),-Grav_grad(:,2), ...
               Grav_grad(:,3)-ZDEM_griddedInterpolant(Grav_grad(:,1),Grav_grad(:,2))-0.5)-...
               GGM_Gravity_griddedInterpolant(Grav_grad(:,1),-Grav_grad(:,2), ...
               Grav_grad(:,3)-ZDEM_griddedInterpolant(Grav_grad(:,1),Grav_grad(:,2))+0.5));
    
    Grav_grad(:,4) = Grav_grad(:,4) - GGM_GravityGradient; 
    
    %plotCustomScatter(Grav_grad(:,1),Grav_grad(:,2),GGM_GravityGradient,GRID_PARA,2*GRID_PARA.buffer,'GGMgravityGradient','mGal/m',OUTPUT_PARA.plotsFolder)
    
    %plotCustomScatter(Grav_grad(:,1),Grav_grad(:,2),Grav_grad(:,4),GRID_PARA,2*GRID_PARA.buffer,'LSCinputGravityGradient','mGal/m',OUTPUT_PARA.plotsFolder)
end
% Run least squares collocation in blocks
block_counter = 0 ; % Increment the counter

for LATsi=GRID_PARA.MAXLAT:-GRID_PARA.STEP:GRID_PARA.MINLAT
for LONGsi=GRID_PARA.MINLONG:GRID_PARA.STEP:GRID_PARA.MAXLONG
    block_counter = block_counter + 1; % Increment the counter
    disp('...................................................................')
    disp(['Working on block ',num2str(block_counter),' centred at: ',num2str(LONGsi),', ',num2str(LATsi)])
    if block_counter >= 2
        break;  % This will exit the loop when the counter reaches the limit
    end
    % Set nodes boundary
    Nodes.Latmin=LATsi;
    Nodes.Latmax=LATsi;
    Nodes.Longmin=LONGsi;
    Nodes.Longmax=LONGsi;

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
    Gravdatout.Long=Gravdatout_m(:,1);
    Gravdatout.Lat=Gravdatout_m(:,2);
    Gravdatout.ortho_H=Gravdatout_m(:,3);
    Gravdatout.Grav_Anom=Gravdatout_m(:,4);
    Gravdatout.STD=Gravdatout_m(:,5);
    %Gravdatout.type=Gravdatout_m(:,6);

    GravGraddatout_m=[];
    if GRAV_GRAD_PARA.avail
        disp('Extracting gravity gradient residual signal in block')
        INOUT_grav_grad=inpolygon(Grav_grad(:,1),Grav_grad(:,2),CoordsMM(:,1),CoordsMM(:,2));
        GravGraddatout_m=Grav_grad(INOUT_grav_grad==1,:);% Format/extracts the airborne gravity data 
        % Assign to vairables to be used thorugh the rest of the code
        GravGraddatout.Long=GravGraddatout_m(:,1);
        GravGraddatout.Lat=GravGraddatout_m(:,2);
        GravGraddatout.ortho_H=GravGraddatout_m(:,3);
        GravGraddatout.Grav_Grad_Anom=GravGraddatout_m(:,4);
        GravGraddatout.STD=GravGraddatout_m(:,5);
    end
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
        
    disp('Pre-computing RTM gradiometry covariance function') 
    CCov_dgdg_int_fun_RTM = precomputeCovarianceFunction('cov_dgdg',RadiusBjerhammar,COV_PARA_RTM.width,COV_PARA_RTM.res,COV_PARA_RTM.A,COV_PARA_RTM.B,COV_PARA_RTM.N,COV_PARA_RTM.M);
    CCov_dgg_int_fun_RTM = precomputeCovarianceFunction('cov_dgg',RadiusBjerhammar,COV_PARA_RTM.width,COV_PARA_RTM.res,COV_PARA_RTM.A,COV_PARA_RTM.B,COV_PARA_RTM.N,COV_PARA_RTM.M);
    CCov_dgt_int_fun_RTM = precomputeCovarianceFunction('cov_dgt',RadiusBjerhammar,COV_PARA_RTM.width,COV_PARA_RTM.res,COV_PARA_RTM.A,COV_PARA_RTM.B,COV_PARA_RTM.N,COV_PARA_RTM.M);
    
    % Compute normal gravity at the computation points on the tile

    NormalGravity=AbsoluteGravityEquator_mgal*(1+NormalGravityConstant*(sin(deg2rad(GRID_REF_dat(:,2))).^2) ...
        )./sqrt(1-EarthEccentricitySquared*(sin(deg2rad(GRID_REF_dat(:,2))).^2));
   
    disp(['The size of the LSC matrix is going to be: ',num2str(length(Gravdatout.Long)),...
        ' x ',num2str(length(Gravdatout.Long))])
   
    if ~isempty(Gravdatout.Long)
        disp('... Interpolating gravimetry covariance matrices')
       
        % Auto-covariance of gravity anomaly at DEM points
        ACOVggRTM_dem = interpolateCovarianceFunction(...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), ...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), CCov_gg_int_fun_RTM);
    
%         % Cross-covariance of gravity anomaly at gravity and DEM points
%         CCOVggRTM_dem_grav = interpolateCovarianceFunction(...
%             GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
%             RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), ...
%             Gravdatout.Long, Gravdatout.Lat, RadiusBjerhammar + Gravdatout.ortho_H, CCov_gg_int_fun_RTM);
%     
%         % Cross-covariance of gravity anomaly and potential at gravity and DEM points
%         CCOVgtRTM_dem_grav = interpolateCovarianceFunction(...
%             GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
%             RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), ...
%             Gravdatout.Long, Gravdatout.Lat, RadiusBjerhammar + Gravdatout.ortho_H, CCov_gt_int_fun_RTM);
        
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

    if ~isempty(GravGraddatout_m)
        
        disp('... Interpolating gradiometry covariance matrices')
        % Auto-covariance of gravity gradient at gravity gradient points
        ACOV_dgdg_grad = interpolateCovarianceFunction(...
                GravGraddatout.Long, GravGraddatout.Lat, RadiusBjerhammar + GravGraddatout.ortho_H,...
                GravGraddatout.Long, GravGraddatout.Lat, RadiusBjerhammar + GravGraddatout.ortho_H, CCov_dgdg_int_fun_RTM);
       
        % Cross-covariance of gravity gradient and potential at gravity gradient and DEM points
        CCOV_dgt_grad_dem = interpolateCovarianceFunction(...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2), ...
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1), GRID_REF_dat(:,2)), ...
            GravGraddatout.Long, GravGraddatout.Lat, RadiusBjerhammar + GravGraddatout.ortho_H, CCov_dgt_int_fun_RTM);
    
        % Cross-covariance of gravity gradient and gravity at DEM points  
        CCOV_dgg_grad_dem=interpolateCovarianceFunction(...
            GRID_REF_dat(:,1), GRID_REF_dat(:,2),...     
            RadiusBjerhammar + ZDEM_griddedInterpolant(GRID_REF_dat(:,1),GRID_REF_dat(:,2)),...
            GravGraddatout.Long, GravGraddatout.Lat, RadiusBjerhammar + GravGraddatout.ortho_H,CCov_dgg_int_fun_RTM);
    
            disp('... Solving LSC matrix equations for Gradiometry') 
     
            [LSCResidualGeoid,errorMatrixGeoid,residualFreeAirGravityAnomaly,residualBouguerGravityAnomaly,...
                errorMatrixGravity]=solveGradientLSCmatrix(GRAV_PARA,GRAV_GRAD_PARA,GravGraddatout,GRID_REF_dat,density, ...
                ACOV_dgdg_grad,ACOVggRTM_dem,CCOV_dgg_grad_dem, ...
                CCOV_dgt_grad_dem,ACOVttRTM_dem,ACOVggFay_dem,CCOVgtFay_dem,RTM_Correction_function,LWLBouguer_Slab_Function,ZDEM_griddedInterpolant,NormalGravity);     
 
    end
    disp('Done')    

   % Create weights for grid blending 

    inputWeightVector=INOUT.*INOUTinner;

    filterWeights = createGridWeights(inputWeightVector, DEM_PARA, GRID_PARA);

    % Save the data to a tile 

    Dataset_save.weights=filterWeights;

    Dataset_save.res_geoid=INOUT*0;
    Dataset_save.res_geoid(INOUT==1)=LSCResidualGeoid;

    Dataset_save.pot_error=INOUT*0;
    Dataset_save.pot_error(INOUT==1)=abs(errorMatrixGeoid);
    
    Dataset_save.res_grav=INOUT*0;
    Dataset_save.res_grav(INOUT==1)=residualFreeAirGravityAnomaly;

    Dataset_save.res_grav_Bouguer=INOUT*0;
    Dataset_save.res_grav_Bouguer(INOUT==1)=residualBouguerGravityAnomaly;

    Dataset_save.grav_error=INOUT*0;
    Dataset_save.grav_error(INOUT==1)=abs(errorMatrixGravity);

    Dataset_save.COV_PARA_RTM=COV_PARA_RTM;
    Dataset_save.COV_PARA_Faye=COV_PARA_Faye;
     
    save([OUTPUT_PARA.Tiles_dir_name,'/Tile',num2str(LONGsi),'_',num2str(abs(LATsi)),'.mat'],'Dataset_save')
    
    else
    disp(('No data in block'))
    end
end
end
%
end