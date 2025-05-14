function [G1, G1_p1, G1_p2] = fftMolodenskyG1(H, FA, Latm, Latmi, Longm, res, sphericalCap)
% FFTG1Deg computes the G1 topographic correction term using FFT-based convolution.
%
% INPUTS:
%   H           : Digital Elevation Model (DEM) in meters
%   FA          : Gravity anomaly values (same size as H), in mGal
%   Latm, Latmi : Latitude grids (same size as H), in degrees
%   Longm       : Longitude grid (same size as H), in degrees
%   res         : Grid resolution in degrees
%   sphericalCap: Radius of spherical cap in degrees
%
% OUTPUT:
%   G1          : Computed G1 correction term in mGal
%
% Note: Missing data should be represented as -9999 in both H and FA
constants
% Preprocessing
disp('Preprocessing input data...')

H(H == -9999) = 0;
FA(FA == -9999) = 0;
FA = FA * 1e-5; % Convert mGal to SI units

% Compute spherical distance (Haversine formula)
dLat = (Latmi - Latm) * pi / 180;
dLong = (mean(Longm(:)) - Longm) * pi / 180;
sinPSIon2 = sqrt(sin(dLat / 2).^2 + cos(Latm * pi / 180) .* cos(Latmi * pi / 180) .* sin(dLong / 2).^2);
distanceFromCP = haversine(Latmi, Longm, Latm, Longm);
l = 2 * ae * sinPSIon2; % Chord length
L = 2 * asin(sinPSIon2); % Angular distance in radians

% Compute integration kernel
Delta = res;
l3 = l.^3;

% Kernel weighting function
W = (l3 ./ (2 * Delta)) .* (((l + Delta / 2).^2 - (l - Delta / 2).^2) ./ ((l + Delta / 2).^2 .* (l - Delta / 2).^2));
K = W ./ l3;

% Apply spherical cap limits
K(L > sphericalCap * pi / 180) = 0;
K(L < res * pi / 180) = 0;

% FFT operations
disp('Computing FFTs...')
FK   = fft2(fftshift(K));
FFA  = fft2(FA);
FHFA = fft2(H .* FA);

% Compute G1 term
disp('Calculating G1 correction...')
areaElement = ((1e5 * ae^2) / (2 * pi)) *(res * pi / 180)^2;
G1 =    (ifft2(FHFA .* FK) - H .* ifft2(FFA .* FK)) * areaElement;
G1_p1 = (ifft2(FHFA.*FK))* areaElement;
G1_p2 = (-1.*ifft2(FFA.*FK))* areaElement;

end
