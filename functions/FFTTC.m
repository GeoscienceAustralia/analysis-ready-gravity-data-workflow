function [TC]=FFTTC(H,Latm,Longm,res,sphericalCap)
% Import the DEM esri ascii data file
%%
tic
H(H==-9999)=0;
H(H<0)=0;
% Extract submatrxices of the DEM heights and location grids 
% and convert the locations to metres
Xm=Longm*111319.9*cos(mean(mean(Latm))*pi/180);
Ym=Latm*111319.9;
% Set up the integration kernel l3
Mlat=mean(mean(Latm));
%%
%% Check if data is all zeros if so jump to the end
if sum(sum(H==0))==numel(H)
   TC=zeros(size(H));
else
%% Calculate integration kernel    
disp('Setting up integration kernel')
l=sqrt((Xm-mean(mean(Xm))).^2+(Ym-mean(mean(Ym))).^2);
clear Xm Ym Longm Latm
Delta=111319.9*res;
l3=l.^(3);
% kernel weighting 
W=(l3./(2*Delta)).*(((l+Delta/2).^2-(l-Delta/2).^2)./(((l+Delta/2).^2).*((l-Delta/2).^2)));
K=W./l3;
% apply a spherical cap to zero the kernel off outside some set radius
K(l>sphericalCap*111319.9)=0;
K(l<111319.9*res)=0;
% Free up some memory
clear l l3 W dz
%% compute 2D fft for the Heights, Heights^2 and kernel
disp('FFT kernel for TC')
FK=fft2(fftshift(K));
clear K
disp('FFT H^2')
FH2=fft2(H.^2);
disp('FFT H')
FH=fft2(H);
disp('FFT ones')
FOne=fft2(ones(size(FK)));
%% Compute the terrain correction in mGal for the whole strip
disp('Computing TC')
TC=abs((((10^5)*2.67*6.6720e-08)/2)*(ifft2(FK.*FH2)...
    -2*H.*ifft2(FK.*FH)+(H.^2).*ifft2(FK.*FOne))...
    *((111319.9*res)*(111319.9*cos(Mlat*pi/180)*res)));
disp(['time taken: ',num2str(toc)])
%% do some plotting and save the output
disp('Done')
end