% This file contains some basic constants:
%  - physical
%  - geodetic (GRS80)
%
% Nico Sneeuw
% Last updated by Neda Darbeheshti
% Geoscience Australia, 2023-11.

% % Some numbers (1986 recommended values, see Physics Today, 08/94, part 2)
% clight = 2.99792458e8;	    % speed of light [m/s]
% G      = 6.67259e-11;         % gravitational constant [m^3 /(kg s^2)]
% au     = 149.597870691e9;     % astronomical unit [m]
% 
% % GRS80 defining constants:
% ae     = 6378136.3;           % semi-major axis of ellipsoid [m]
% GM     = 3.986004415e14;      % geocentric grav. constant [m^3 / s^2]
% J2     = 1.08263e-3;          % earth's dyn. form factor (= -C20 unnormalized)
% Omega  = 7.292115e-5;         % mean ang. velocity [rad/s]
% 
% % GRS80 derived constants:
flattening   = 1/298.257222101; % flattening      
% J4     = -0.237091222e-5;     % -C40 unnormalized
% J6     =  0.608347e-8;        % -C60 unnormalized
% J8     = -0.1427e-10;         % -C80 unnormalized

% mGal   = 10^(-5) m / s^2
% 1 E    = 10^(-9) m / s^2 / m

% % GRS80 constants 

EarthMajorAxis = 6378.137;                                % Earth Major Axis [km]
%computeLSC
EarthMinorAxis = 6356.752;                                % Earth Minor Axis [km]
%computeLSC
EarthRadius = 6371000;                                    % Earth Radius [m]
EarthEccentricitySquared = 0.00669438002290;                % Earth Eccentricity Squared
%computeLSC

% % physical constants

NormalGravityConstant = 0.001931851353;                     % Normal Gravity Constant
AbsoluteGravityEquator_mgal = 9.7803267715*(10^5);          % Absolute Gravity Equator [mgal]
AbsoluteGravityPole_mgal =    9.8321849378*(10^5);          % Absolute Gravity Pole [mgal]
GravToCentrifugalRatio_Equator = 0.00344978600308;          % the ratio between the gravitational and centrifugal forces
                                                            % defined as (w^2*a^2*b)/GM
%computeLSC
BouguerConstant=0.0419;                                   % Bouguer Constant
%computeLSC,computeTerrainEffect
bigG = 6.6720e-08;                                        % gravitational constant [dyn⋅cm2⋅g−2]
%computePrismGravity 

% % calculation constants

deg2meter=111319.9;                                        % convert from degree to meter 
%computeTerrainCorrection