function [G1, G1_p1, G1_p2] = fftMolodenskyG1(DEM, freeAirGravity, latDEM, longDEM, res, sphericalCap)
% FFTG1Deg computes the G1 topographic correction term using FFT-based convolution.
%
% INPUTS:
%   DEM            : Digital Elevation Model (DEM) in meters
%   freeAirGravity : Gravity anomaly values (same size as H), in mGal
%   Latm, Latmi    : Latitude grids (same size as H), in degrees
%   Longm          : Longitude grid (same size as H), in degrees
%   res            : Grid resolution in degrees
%   sphericalCap   : Radius of spherical cap in degrees
%
% OUTPUT:
%   G1          : Computed G1 correction term in mGal
%
% Note: Missing data should be represented as -9999 in both H and FA
% Last updated by Neda Darbeheshti
% Geoscience Australia, 2025-05.
constants
% Preprocessing input data...
DEM(DEM == -9999) = 0;
freeAirGravity(freeAirGravity == -9999) = 0;
freeAirGravity = freeAirGravity * 1e-5; % Convert mGal to SI units
% Compute spherical distance (Haversine formula)
Latmi = mean(latDEM(:,1));
LatmiRad = deg2rad (Latmi);
LatmRad = deg2rad (latDEM);
LongmRad = deg2rad (longDEM);
L = haversine(LatmiRad, mean(LongmRad(:)), LatmRad, LongmRad);
littleL=2*ae*sin(L/2);
% Compute integration kernel
l3 = littleL.^3;

% Kernel weighting function
KernelWeight = (l3 ./ (2 * res)) .* (((littleL + res / 2).^2 - (littleL - res / 2).^2) ./ ((littleL + res / 2).^2 .* (littleL - res / 2).^2));
K = KernelWeight ./ l3;

% Apply spherical cap limits
K(L > deg2rad (sphericalCap)) = 0;
K(L < deg2rad (res)) = 0;

% FFT operations
FK   = fft2(fftshift(K));
FFA  = fft2(freeAirGravity);
FHFA = fft2(DEM .* freeAirGravity);

% Compute G1 term
areaElement = ((1e5 * ae^2) / (2 * pi)) *(deg2rad (res))^2;
G1 =    (ifft2(FHFA .* FK) - DEM .* ifft2(FFA .* FK)) * areaElement;
G1_p1 = (ifft2(FHFA.*FK))* areaElement;
G1_p2 = (-1.*ifft2(FFA.*FK))* areaElement;

end
