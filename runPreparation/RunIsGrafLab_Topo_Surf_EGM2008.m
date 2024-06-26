%% To be used for EGM2008 geoid values and gravity anomalies.
%
% Main parameters to be defined:
%                      - Grid Parameters
% Main functions:
%             - isGrafLabJM: 
% Written by Jack McCubbine
% Last updated by Neda Darbeheshti
% Geoscience Australia, 2024-05.
clear 
close all
%% Add the path to the function files.
addpath('functions');
constants                                       % load constants
%% Set Ellipsoidal parameters
EarthMajorAxis=EarthMajorAxis*1000;
EarthMinorAxis=EarthMinorAxis*1000;
AbsoluteGravityEquator_mgal=AbsoluteGravityEquator_mgal*(10^(-5));
AbsoluteGravityPole_mgal=AbsoluteGravityPole_mgal*(10^(-5));
%% Set some parameters
plotsFolder=['Outputs/plots/',date,'grafLab'];
outputFile='Data/processedData/GOCE_N300_For_Gridded_Int24_2height.mat';%EGM2008_For_Gridded_Int24.mat
GGMfile='Data/GGM/GO_CONS_GCF_2_DIR_R6.gfc';%egm2008zerotide.gfc'
DEMfile='Data/DEM/AUSDEM1m.xyz';
resolution='1/60';
Longmin=93;
Longmax=174;
Latmin=-61;
Latmax=-8;
H_STEP=2000;%250 for originnal 9 layer hight 
H_MAX=2000;
%% Sort out DEM
ncols=round((Longmax-Longmin)*1./str2num(resolution))+1;
nrows=round((Latmax-Latmin)*1./str2num(resolution))+1;
SMT=importdata(DEMfile);
Z=reshape(SMT(:,3),ncols,nrows)';
% Compute and save 3 versions of it
Z(Z<-100)=0;
save('Data/processedData/3sZ.mat','Z');
%% Get ellipsoidal gravity 
% Gravity disturbance sa 21
% Deflection of the vertical xi 3 %north south deflection 
% Disturbing potential 5
% Height anomaly ell 22
% Second radial derivative of disturbing potential 24
isGrafLabJM(0,GGMfile,Latmin,Latmax,Longmin,Longmax,resolution,0,'Data/processedData/3sZ.mat',21,3,5,22,'Data/processedData/GetDg_ell_Topo_Surf','')
% Convert to ellipsoidal gravity anomaly
GL1=importdata('Data/processedData/GetDg_ell_Topo_Surf.mat');
heightAnomaly=GL1(:,7);
Ham=reshape(heightAnomaly,3181,4861);
Z2=Z+Ham(end:-1:1,:);
save('Data/processedData/3sZ.mat','Z2');
%% Get ellipsoidal gravity 
close all
isGrafLabJM(0,GGMfile,Latmin,Latmax,Longmin,Longmax,resolution,0,'Data/processedData/3sZ.mat',21,3,5,22,'Data/processedData/GetDg_ell_Topo_Surf','')
% Convert to ellipsoidal gravity anomaly
GL2=importdata('Data/processedData/GetDg_ell_Topo_Surf.mat');
Hao=GL2(:,7);
%% Gravity anomaly on the ellipsoid + 4000
disp('Calculating on the Surface plus H')
Hcount=0;
Deltagms=zeros(3181,4861, Hcount);
Hams=zeros(3181,4861, Hcount);
Trrs=zeros(3181,4861, Hcount);
for H=0:H_STEP:H_MAX
    Hcount=Hcount+1;
%     if Hcount >= 2
%         break;  % This will exit the loop when the counter reaches the limit
%     end
Ham=reshape(Hao,3181,4861);
Z2=Z+Ham(end:-1:1,:)+H;
save('Data/processedData/3sZ.mat','Z2');
%%
close all
isGrafLabJM(0,GGMfile,Latmin,Latmax,Longmin,Longmax,resolution,0,'Data/processedData/3sZ.mat',21,3,5,22,'Data/processedData/GetDg_ell_Topo_Surf','')
close all
isGrafLabJM(0,GGMfile,Latmin,Latmax,Longmin,Longmax,resolution,0,'Data/processedData/3sZ.mat',21,3,24,22,'Data/processedData/GetDg_ell_Topo_Surf24','')
% Convert to ellipsoidal gravity anomaly
iGL=importdata('Data/processedData/GetDg_ell_Topo_Surf.mat');
lat=iGL(:,1);
lat_rad=deg2rad(lat);
long=iGL(:,2);
height=iGL(:,3);
deltag_sa=iGL(:,4)*10^-5;
verticalDeflection=deg2rad(iGL(:,5)/3600);
potential=iGL(:,6);
heightAnomaly=iGL(:,7);
iGL24=importdata('Data/processedData/GetDg_ell_Topo_Surf24.mat');
Trr=iGL24(:,6);
%% Get geocentric radius
[X,Y,Zr]=geodetic2ecef(lat_rad,0.*zeros(length(lat),1),height.*ones(length(lat),1),[EarthMajorAxis sqrt(EarthEccentricitySquared)]);  
r=sqrt(X.*X+Y.*Y+Zr.*Zr); %Radius
%% Gamma, page 27 McCubbine 2016 thesis
gamma_ell=(EarthMajorAxis*AbsoluteGravityEquator_mgal*(cos(lat_rad).^2)+EarthMinorAxis*AbsoluteGravityPole_mgal*(sin(lat_rad).^2))./sqrt((EarthMajorAxis^2*(cos(lat_rad).^2)+EarthMinorAxis^2*(sin(lat_rad).^2)));
gamma=gamma_ell.*(1-(2*height/EarthMajorAxis).*(1+flattening+GravToCentrifugalRatio_Equator-2*flattening*(sin(lat_rad).^2))+3*(height.^2)./(EarthMajorAxis.^2));
dgamma_ell_dh=-(2*gamma./EarthMajorAxis).*(1+flattening+GravToCentrifugalRatio_Equator-2*flattening*(sin(lat_rad).^2));
%% Partial derivatives
pLat_rad=deg2rad(90+abs(lat));
v=EarthMajorAxis./sqrt(1-EarthEccentricitySquared*cos(pLat_rad).^2);
partialrpartialh=(v.*(1-EarthEccentricitySquared*(cos(pLat_rad).^2))+height)./r;
partialthetapartialh=-EarthEccentricitySquared.*v.*tan(pLat_rad)./((v*(1-EarthEccentricitySquared)+height).^2+((v+height).^2).*(tan(pLat_rad).^2));
%% Get Detlag_ell
Deltag_ell=((partialrpartialh).*deltag_sa-(partialthetapartialh.*r.*gamma.*verticalDeflection)+(1./gamma).*dgamma_ell_dh.*potential)*10^5;

% Convert to matrices and plot
Deltagm=reshape(Deltag_ell,3181,4861);
Deltagm=Deltagm(end:-1:1,:); % This effectively flips the matrix vertically.

Ham=reshape(heightAnomaly,3181,4861);
Ham=Ham(end:-1:1,:);

Trrm=reshape(Trr,3181,4861);
Trrm=Trrm(end:-1:1,:);

Trrs(:,:,Hcount)=Trrm;
Deltagms(:,:,Hcount)=Deltagm;
Hams(:,:,Hcount)=Ham;

Latm=reshape(lat,3181,4861);
Latm=Latm(end:-1:1,:);

Longm=reshape(long,3181,4861);
Longm=Longm(end:-1:1,:);

figure(Hcount)
clf
hold on
imagesc(Longm(1,:),Latm(:,1),Deltagms(:,:,Hcount))
title('Deltagms')
colorbar
drawnow

figure(Hcount+10)
clf
hold on
imagesc(Longm(1,:),Latm(:,1),Trrs(:,:,Hcount))
title('Trrs')
colorbar
drawnow
end
% prepare to save
[y,x,z]=ndgrid(-Latm(:,1),Longm(1,:),0:H_STEP:H_MAX);
P = [2 1 3];
x = permute(x, P);
y = permute(y, P);
z = permute(z, P);
Deltagms = permute(Deltagms, P);
Hams = permute(Hams, P);
Trrs = permute(Trrs, P);

Forin.x=x;
Forin.y=y;
Forin.z=z;
Forin.g=Deltagms;
Forin.zeta=Hams;
Forin.trr=Trrs;

save(outputFile,'-struct','Forin')

% Test it out
Hams_griddedInterpolant=griddedInterpolant(x,y,z,Hams);
Hams_plot=Hams_griddedInterpolant(Longm(:),-Latm(:),1000+0*Latm(:))-Hams_griddedInterpolant(Longm(:),-Latm(:),2000+0*Latm(:));

Trrs_griddedInterpolant=griddedInterpolant(x,y,z,Trrs);
Trrs_plot=Trrs_griddedInterpolant(Longm(:),-Latm(:),1000+0*Latm(:))-Trrs_griddedInterpolant(Longm(:),-Latm(:),2000+0*Latm(:));

Deltagms_griddedInterpolant=griddedInterpolant(x,y,z,Deltagms);
Deltagms_plot=Deltagms_griddedInterpolant(Longm(:),-Latm(:),1000+0*Latm(:))-Deltagms_griddedInterpolant(Longm(:),-Latm(:),2000+0*Latm(:));

% common variables for plotting
Coastline=importdata('Data\COASTLINE\CoastAus.mat');
axisLimits.latMeanCosine=abs(cos(deg2rad(mean([Latmin Latmax]))));
axisLimits.lonMinLimit=Longmin;
axisLimits.lonMaxLimit=Longmax;
axisLimits.latMinLimit=Latmin;
axisLimits.latMaxLimit=Latmax;

% plots
figure('Name','Height','NumberTitle','off');
clf
hold on
imagesc(Longm(1,:),Latm(:,1),reshape(Hams_plot,3181,4861))
customizeMap('Height','m',Coastline,axisLimits)

figure('Name','Trr','NumberTitle','off');
clf
hold on
imagesc(Longm(1,:),Latm(:,1),reshape(Trrs_plot,3181,4861))
customizeMap('Trr','mGal/m',Coastline,axisLimits)

figure('Name','Delta g','NumberTitle','off');
clf
hold on
imagesc(Longm(1,:),Latm(:,1),reshape(Deltagms_plot,3181,4861))
customizeMap('Delta g','mGal',Coastline,axisLimits)


