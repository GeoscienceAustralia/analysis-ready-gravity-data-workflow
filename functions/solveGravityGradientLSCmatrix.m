function [res_Pot_g,Geoid_Error,residualFreeAirGravityAnomaly,res_grav_g_Bouguer,Grav_Error]= solveGravityGradientLSCmatrix(GRAV_PARA,GRAV_GRAD_PARA,Gravdatout,GravGraddatout,GRID_REF_dat,density, ...
            ACOV_gg_grav,ACOV_dgdg_grad,ACOV_dgg_grad_grav,ACOV_gg_dem,CCOV_gg_grav_dem,CCOV_dgg_grad_dem,CCOV_gt_grav_dem, ...
            CCOV_dgt_grad_dem,ACOV_tt_dem,ACOV_gg_dem_Faye,CCOV_gt_dem_Faye,RTM_Correction_function,LWLBouguer_Slab_Function,ZDEM_Function,NormalGravity)

    % performs LSC twice for combination of gravity and gravity gradient
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

    LSCMAT_comb=[ACOV_dgdg_grad+diag(GravGraddatout.STD.^2+GRAV_GRAD_PARA.TypeB^2),ACOV_dgg_grad_grav';...
            ACOV_dgg_grad_grav, ACOV_gg_grav+diag(Gravdatout.STD.^2+GRAV_PARA.TypeB^2)];

    % inverse of covariance matrix
    inverseCovarianceMatrix=LSCMAT_comb\eye(size(LSCMAT_comb));
    
    %middleMatrix=inverseCovarianceMatrix*[GravGraddatout.Grav_Grad_Anom;Gravdatout.Grav_Anom];
    
    % LSC residual gravity anomly
    %res_grav_g=[CCOV_dgg_grad_dem,CCOV_gg_grav_dem]*middleMatrix;

    res_grav_g=[CCOV_dgg_grad_dem,CCOV_gg_grav_dem]*(LSCMAT_comb\[GravGraddatout.Grav_Grad_Anom;Gravdatout.Grav_Anom]);
     

    % LSC residual gravity anomly - RTM correction
    residualFreeAirGravityAnomaly=res_grav_g-RTM_Correction_function(GRID_REF_dat(:,1),GRID_REF_dat(:,2));

    % LSC residual gravity anomly - long wavelength Bouguer 
    res_grav_g_Bouguer=res_grav_g-LWLBouguer_Slab_Function(GRID_REF_dat(:,1),GRID_REF_dat(:,2));

    % Error Covariance Matrix
    Grav_Error=sqrt(diag(ACOV_gg_dem-[CCOV_dgg_grad_dem,CCOV_gg_grav_dem]*inverseCovarianceMatrix*[CCOV_dgg_grad_dem,CCOV_gg_grav_dem]')); 
    
    % residualBouguerGravityAnomaly + Bouguer formula eq.(2.14) Mccubbine(2016)
    Faye=res_grav_g_Bouguer+ZDEM_Function(GRID_REF_dat(:,1),GRID_REF_dat(:,2))*BouguerConstant*density;% This is hard coded density, revisit.

    % LSC residual geoid
    res_Pot_g=(CCOV_gt_dem_Faye*((ACOV_gg_dem_Faye+diag(Grav_Error.^2+GRAV_PARA.Grav_Faye_TypeB^2))\Faye))./NormalGravity; 
   
    % error matrix
    Geoid_Error=sqrt(diag(ACOV_tt_dem-[CCOV_dgt_grad_dem,CCOV_gt_grav_dem]*inverseCovarianceMatrix*[CCOV_dgt_grad_dem,CCOV_gt_grav_dem]'))./NormalGravity;
end