
constants
% Covariance function parameters
COV_PARA.Compute_Empircal_COV_Dec=3; % Decimation factor for empirical covariance estimates. e.g. 1 is no decimation, 2 drops 50% of the data etc. see sph_empcov for logic.
COV_PARA.Fit_Empircal_COV='auto';%'auto';% process to fit covariance N & M function values 'man' for manual to fit them on the cmd line or 'auto' , '' to just use what you supply here.
COV_PARA.FitEmpiricalCOVNSearch=[21600,1,21600];%10800 is not working:speckelly pattern% Start, step, stop parameter sweep values for N parameter - if auto
COV_PARA.FitEmpiricalCOVMSearch=[200,20,300];% Start, step, stop parameter sweep values for M parameter - if auto
COV_PARA.N=10800;% max Legender polynonial of cov func. 
COV_PARA.M=200;% min Legender polynonial of cov func. 
COV_PARA.width=3;% Size of precomputed cov function in degrees - must be larger the the distance between any two points on a tile. 
COV_PARA.res=30/3600; % Resolution of the covariance function
COV_PARA.COV_COMPUTED_Tilewise=true;% This recomputes the covariance function for each tile.
COV_PARA.Airbornedataonly=false;%Only use airborne data in establishing Covariance parameters - good to use if we are using EGM2008 as the references as terrestrial data are not independent.
COV_PARA.COVPlot=false;% true plots progress, false turns this off.

% radius of the Bjerhammar sphere
phi=deg2rad(mean([-31.5,-33.5]));
RadiusBjerhammar= EarthMajorAxis*EarthMinorAxis/sqrt((EarthMajorAxis*sin(phi)).^2+(EarthMinorAxis*cos(phi)).^2)*10^3;% Pajama sphere radius.

coefficientA = 84.2143;
coefficientB = -294.2557;

 CCov_tt_int_fun_RTM = precomputeCovarianceFunctionTest4all('cov_tt', RadiusBjerhammar, ...
    COV_PARA.width, COV_PARA.res, coefficientA, coefficientB, COV_PARA.N, COV_PARA.M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COVA([], 1)