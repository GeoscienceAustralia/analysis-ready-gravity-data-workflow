%%%%%%%%%%%%%%%%%%%%%%%% new way to make the covariance matrix
constants                                       % load constants
phi=deg2rad(mean(GRID_REF(:,2)));
RadiusBjerhammar= EarthMajorAxis*EarthMinorAxis/sqrt((EarthMajorAxis*sin(phi)).^2+(EarthMinorAxis*cos(phi)).^2)*10^3;% Pajama sphere radius.

CCov_tt_int_fun_RTM=precomputeCovarianceFunction('cov_tt',RadiusBjerhammar,COV_PARA.width,COV_PARA.res,sigma2,bestFitCoeff,COV_PARA.N,COV_PARA.M);

% Auto-covariance of potential at DEM points
ACOVttRTM_lev = interpolateCovarianceFunction(...
Lev(:,1), Lev(:,2), ...
RadiusBjerhammar + ZDEM_griddedInterpolant(Lev(:,1), Lev(:,2)), ...
Lev(:,1), Lev(:,2), ...
RadiusBjerhammar + ZDEM_griddedInterpolant(Lev(:,1), Lev(:,2)), CCov_tt_int_fun_RTM,OUTPUT_PARA,'ACOVttGPSlevelling m^4/s^4',1);
figure(3)
imagesc(ACOVttRTM_lev)