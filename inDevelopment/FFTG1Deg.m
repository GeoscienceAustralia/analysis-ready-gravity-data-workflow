function [G1]=FFTG1Deg(H,FA,Latm,Latmi,Longm,res,sphericalCap)
%FFTG1(filename,filenameFA,minlat,maxlat,minlong,maxlong,res,sphericalCap,outfile,PadS,PadN,PadW,PadE)
% Import the DEM esri ascii data file
%%

H(H==-9999)=0;
FA(FA==-9999)=0;
FA(isnan(FA))=0;
H(isnan(H))=0;
FA=FA*(10^-5);
%% Calculate integration kernel    
disp('Setting up integration kernel')
R=0.6378136460E+07;
sinPSIon2=sqrt(sin((Latmi*pi/180-Latm*pi/180)/2).^2+(sin((mean(mean(Longm*pi/180))-Longm*pi/180)/2).^2).*cos(Latm*pi/180)*cos(Latmi*pi/180));
size(sinPSIon2)
l=2*R*sinPSIon2;
L=2*asin(sinPSIon2);
clear Longm Latm
Delta=res;
l3=l.^(3);
% kernel weighting 
W=(l3./(2*Delta)).*(((l+Delta/2).^2-(l-Delta/2).^2)./(((l+Delta/2).^2).*((l-Delta/2).^2)));
K=W./l3;
% apply a spherical cap to zero the kernel off outside some set radius
K(L>sphericalCap*pi/180)=0;
K(L<res*pi/180)=0;
% Free up some memory
clear l l3 W dz
%% compute 2D fft
disp('FFT kernel for G1')
FK=fft2(fftshift(K));
clear K
disp('FFT FA')
FFA=fft2(FA);
disp('FFT H.*FA')
FHFA=fft2(H.*FA);
%% Compute the G1 term in mGal for the whole strip
disp('Computing G1')
G1=((10^5)*(R^2)/(2*pi))*(ifft2(FHFA.*FK)-H.*ifft2(FFA.*FK)).*((res*pi/180)*(res*pi/180));
end







