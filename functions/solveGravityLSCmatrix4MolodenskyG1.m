function [LSCResidualGravity]=solveGravityLSCmatrix4MolodenskyG1(GRAV_PARA,Gravdatout, ...
            ACOVggRTM_grav,CCOVggRTM_dem_grav)
                
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
  
end