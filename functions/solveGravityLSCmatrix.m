function [LSCResidualGeoid,errorMatrixGeoid,residualFreeAirGravityAnomaly,residualBouguerGravityAnomaly,...
            errorMatrixGravity]=solveGravityLSCmatrix(GRAV_PARA,Gravdatout,GRID_REF_dat,density, ...
            ACOVggRTM_grav,ACOVggRTM_dem,CCOVggRTM_dem_grav,CCOVgtRTM_dem_grav,ACOVttRTM_dem,ACOVggFay_dem, ...
            CCOVgtFay_dem,RTM_Correction_function,LWLBouguer_Slab_Function,ZDEM_griddedInterpolant,NormalGravity)
                
    % performs LSC twice for gravity 
    %
    % Input:  GRID_PARA = struct with buffer,buffer2,STEP,MINLONG,MAXLONG,MINLAT,MAXLAT 
    %         Topo_PARA = struct with Corr, TopoPlot,Density,Depth,Rad,RTM[,,] 
    %
    % Output:   LSCResidualGeoid =
    %           errorMatrixGeoid =
    %           residualFreeAirGravityAnomaly =
    %           residualBouguerGravityAnomaly =
    %           errorMatrixGravity =
    % Example: see computeGravimetryGradiometryLSC
    %
    % Main functions
    % - 
    % -  
    % Other functions
    % -  
    % -
    %
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2024-03.
    constants                                       % load constants

    % inverse of auto covariance matrix
    inverseCovarianceMatrix=(ACOVggRTM_grav+diag(Gravdatout.STD.^2+GRAV_PARA.TypeB^2))\eye(size(ACOVggRTM_grav));
    
    % LSC residual gravity anomly
    LSCResidualGravity=CCOVggRTM_dem_grav*inverseCovarianceMatrix*(Gravdatout.Grav_Anom);

    % LSC residual gravity anomly - RTM correction
    residualFreeAirGravityAnomaly=LSCResidualGravity-RTM_Correction_function(GRID_REF_dat(:,1),GRID_REF_dat(:,2));
    
    % LSC residual gravity anomly - long wavelength Bouguer 
    residualBouguerGravityAnomaly=LSCResidualGravity-LWLBouguer_Slab_Function(GRID_REF_dat(:,1),GRID_REF_dat(:,2));
   
    % residualBouguerGravityAnomaly + Bouguer formula eq.(2.14) Mccubbine(2016)
    residualFayeGravityAnomaly=residualBouguerGravityAnomaly+ZDEM_griddedInterpolant(GRID_REF_dat(:,1),GRID_REF_dat(:,2))*BouguerConstant*density;
    
    % Error Covariance Matrix, eq.(7.8) Mccubbine(2016)
    errorMatrixGravity=sqrt(diag(ACOVggRTM_dem-CCOVggRTM_dem_grav*inverseCovarianceMatrix*CCOVggRTM_dem_grav')); 
    
    % LSC residual geoid
    LSCResidualGeoid=(CCOVgtFay_dem*((ACOVggFay_dem + diag(errorMatrixGravity.^2+GRAV_PARA.Grav_Faye_TypeB^2))\residualFayeGravityAnomaly))./NormalGravity; 
   
    % error matrix, eq.(7.8) Mccubbine(2016)
    errorMatrixGeoid=sqrt(diag(ACOVttRTM_dem - CCOVgtRTM_dem_grav*inverseCovarianceMatrix*CCOVgtRTM_dem_grav'))./NormalGravity;
end