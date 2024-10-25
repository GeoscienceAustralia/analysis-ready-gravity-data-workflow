function fac=computeSecondOrderFreeAirCorrection(lat,ht)
%2nd order Free-Air Correction for WGS-84 Ellipsoid
%reference: GRAV-D General Airborne Gravity Data User Manual, page 23
%Data Inputs
%1.Geodetic latitude: lat
%2.ELLIPSOIDAL height: ht
%Define Ellipsoidal Parameters
a = 6378137; %semi-major axis, WGS-84
g_e = 978032.53359; %equatorial normal gravity, WGS-84
f = 0.00335281066474; %flattening, WGS84
m = 0.00344978650684; %defined as (w^2*a^2*b)/GM, WGS-84

%use GEOCENTRIC latitude, which is defined as: phi
phi = atan(tand(lat).*((1-f).^2));
sphi2 = sin(phi).^2;

c1 = ((2.*g_e)./a);
term1 = (1+f+m-(2.*f.*sphi2)).*ht;
c2 = (3.*g_e./(a.^2));
term2 =ht.^2;

fac = (c1.*term1)-(c2.*term2);


constants                                       % load constants

EarthMajorAxis= EarthMajorAxis*10^3;            % km to m

c1 = ((2.*NormalGravity)./EarthMajorAxis);

term1 = (1+flattening+GravToCentrifugalRatio_Equator-2*flattening*sin(data(:,2)*pi/180).^2).*ht;

fac=(2*(NormalGravity/a).*(1+f+GravToCentrifugalRatio_Equator-2*f*sin(data(:,2)*pi/180).^2).*data(:,3) - (data(:,3).^2)*3.*NormalGravity/a^2);




