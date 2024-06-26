function [res_Pot_g,errorMatrixGeoid,residualFreeAirGravityAnomaly,res_grav_g_Bouguer,...
            errorMatrixGravity]= solveGradientLSCmatrix(GRAV_PARA,GRAV_GRAD_PARA,GravGraddatout,GRID_REF_dat,density, ...
            ACOV_dgdg_grad,ACOV_gg_dem,CCOV_dgg_grad_dem,CCOV_dgt_grad_dem,ACOV_tt_dem,ACOV_gg_dem_Faye, ...
            CCOV_gt_dem_Faye,RTM_Correction_function,LWLBouguer_Slab_Function,ZDEM_Function,NormalGravity)

    % performs LSC twice for gravity gradient
    %
    % Input:  GRID_PARA = struct with buffer,buffer2,STEP,MINLONG,MAXLONG,MINLAT,MAXLAT 
    %         Topo_PARA = struct with Corr, TopoPlot,Density,Depth,Rad,RTM[,,] 
    %
    % Output:   LSCResidualGeoid =
    %           errorMatrixGeoid =
    %           residualFreeAirGravityAnomaly =
    %           residualBouguerGravityAnomaly =
    %           errorMatrixGravity =
    %
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
    inverseCovarianceMatrix=(ACOV_dgdg_grad+diag(GravGraddatout.STD.^2+GRAV_GRAD_PARA.TypeB^2))\eye(size(ACOV_dgdg_grad));

    % LSC residual gravity anomly
    LSCResidualGravity=CCOV_dgg_grad_dem*inverseCovarianceMatrix*(GravGraddatout.Grav_Grad_Anom);
    
    % LSC residual gravity anomly - RTM correction
    residualFreeAirGravityAnomaly=LSCResidualGravity-RTM_Correction_function(GRID_REF_dat(:,1),GRID_REF_dat(:,2));

    % LSC residual gravity anomly - long wavelength Bouguer 
    res_grav_g_Bouguer=LSCResidualGravity-LWLBouguer_Slab_Function(GRID_REF_dat(:,1),GRID_REF_dat(:,2));

    % Error Covariance Matrix
    errorMatrixGravity=sqrt(diag(ACOV_gg_dem-CCOV_dgg_grad_dem*inverseCovarianceMatrix*CCOV_dgg_grad_dem')); 
    
    % residualBouguerGravityAnomaly + Bouguer formula eq.(2.14) Mccubbine(2016)
    Faye=res_grav_g_Bouguer+ZDEM_Function(GRID_REF_dat(:,1),GRID_REF_dat(:,2))*BouguerConstant*density;% This is hard coded density, revisit.

    % LSC residual geoid
    res_Pot_g=(CCOV_gt_dem_Faye*((ACOV_gg_dem_Faye+diag(errorMatrixGravity.^2+GRAV_PARA.Grav_Faye_TypeB^2))\Faye))./NormalGravity; 
   
    % error matrix
    errorMatrixGeoid=sqrt(diag(ACOV_tt_dem - CCOV_dgt_grad_dem*inverseCovarianceMatrix*CCOV_dgt_grad_dem'))./NormalGravity;

    end