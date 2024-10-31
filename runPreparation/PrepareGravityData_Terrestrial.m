close all
clear 
%% Physcial parameters
% a=6378137.d0;				
% a2=a*a;
% b=6356752.3141d0;			
% b2=b*b;		
% f=1.d0/298.257222101d0;
% m=0.00344978600308d0;   
% e2=0.00669438002290d0;
% e=0.00669438002290;
% G=6.67428d-11;
constants                                       % load constants
%% Import the Sandwell and smith gravity anomaly
llf=importdata('Data/GRAVITY/ALTIMETRY/sand311ausgrav.llf');
lle=importdata('Data/GRAVITY/ALTIMETRY/sand311ausgrav.lle');
FreeAirmo=reshape(llf(:,3),4860,3180)';
FreeAir_Errormo=reshape(lle(:,3),4860,3180)';
Latmo=reshape(llf(:,2),4860,3180)';
Longmo=reshape(llf(:,1),4860,3180)';
%% Import DEM
DEM=importdata('Data/DEM/AUSDEM1min.xyz');
Latm_DEM=reshape(DEM(:,2),length(unique(DEM(:,1))),length(unique(DEM(:,2))))';
Longm_DEM=reshape(DEM(:,1),length(unique(DEM(:,1))),length(unique(DEM(:,2))))';
Hm_DEM=reshape(DEM(:,3),length(unique(DEM(:,1))),length(unique(DEM(:,2))))';

% Resample
Latm=Latm_DEM;
Longm=Longm_DEM;
FreeAirm=griddata(Latmo,Longmo,FreeAirmo,Latm,Longm);
% add 1 or 2 to the error
FreeAir_Errorm=griddata(Latmo,Longmo,FreeAir_Errormo,Latm,Longm);
%% Mask out altimetry data onshore (i.e. H>0)
MASK=0*Latm;
H_Function=griddedInterpolant(Longm_DEM',Latm_DEM(end:-1:1,:)',Hm_DEM(end:-1:1,:)');
MASK(H_Function(Longm,Latm)~=0)=1;    
%% Import the raw GADDS gravity data
GADDS=fopen('Data/GRAVITY/GADDS/gravity_may2016.dat');
TEMP=strsplit(fgetl(GADDS));

t=0;
%tic
data=zeros(1800000,9);
for k=1:length(data(:,1))
    try
                   %[Longitude, Latitude, Height, height, N, Gravity, Horizontal Uncertainty, Vertical Uncertainty]
        data(k,:)=[str2double(TEMP{5}),str2double(TEMP{6}),str2double(TEMP{7}),str2double(TEMP{8}),str2double(TEMP{9}),str2double(TEMP{10})/10,str2double(TEMP{19}),str2double(TEMP{20}),str2double(TEMP{23})/10];
        TEMP=strsplit(fgetl(GADDS));
        t=k;
        if mod(k,1000)==0
            toc
            tic
            disp(k)
        end
    catch
    end
end
data(t+1:end,:)=[];
%toc
fclose(GADDS);
% Convert degrees to radians
phi=deg2rad(data(:,2));

% Compute the corrections

NormalGravity = computeNormalGravity (phi);

SecondOrderFreeAirCorrection = computeSecondOrderFreeAirCorrection (phi,data(:,3),NormalGravity);
atmospheric_corr=computeAtmosphericGravityCorrection (data(:,3));

atmospheric_corr1=0.871-1.0298*(10^-4)*data(:,3)+5.3105*(10^-9)*(data(:,3).^2)-2.1642*(10^-13)*(data(:,3).^3);
% Anomaly
FREE_AIR_Anomaly=data(:,6)-NormalGravity+SecondOrderFreeAirCorrection+atmospheric_corr;
FREE_AIR_Anomaly_Error=sqrt(data(:,9).^2+(data(:,7).^2).*(0.000812*sin(2*data(:,2)*pi/180)).^2+((0.3086-0.1119).^2).*((data(:,8).^2)));
disp('Combining datasets and blockmean')
%% Combine the datasets
% Define the flag column
flagColumnTerresterial = 1 * ones(size(data, 1), 1);
flagColumnAltimetry = 2 * ones(size(Longm(MASK==0), 1), 1);
% Add the flag column to the Combined_Grav matrix
Combined_Grav=[Longm(MASK==0),Latm(MASK==0),FreeAirm(MASK==0),flagColumnAltimetry;data(:,1),data(:,2),FREE_AIR_Anomaly-0.0419*2.67*data(:,3),flagColumnTerresterial];
Combined_Error=[Longm(MASK==0),Latm(MASK==0),FreeAir_Errorm(MASK==0),flagColumnAltimetry;data(:,1),data(:,2),FREE_AIR_Anomaly_Error,flagColumnTerresterial];
%% Block average Bougeur
Long1am=round(Combined_Grav(:,1)*60)/60;
Lat1am=round(Combined_Grav(:,2)*60)/60;

flags=Combined_Grav(:,4);
[valsf, uidf, idxf]=unique([Long1am,Lat1am,flags],'rows');

[vals, uid, idx]=unique([Long1am,Lat1am],'rows');

Blockmean_flags = accumarray(idx,Combined_Grav(:,4),[],@mean);
Blockmedian_flags = accumarray(idx,Combined_Grav(:,4),[],@median);
Blockmin_flags = accumarray(idx,Combined_Grav(:,4),[],@min);
Blockmax_flags = accumarray(idx,Combined_Grav(:,4),[],@max);

Blockmean_anom = accumarray(idx,Combined_Grav(:,3),[],@mean);
Blockmean_anom_Err = accumarray(idx,Combined_Error(:,3),[],@mean);
Blockmean_llf=[vals,Blockmean_anom];
Blockmean_Err=abs(sqrt(Blockmean_anom_Err.^2));
disp('done')
%% Create the output file
% Data structure:
% Grav_out=[Longitude,Latitude,Ortho-H,gravity anomaly,uncertainty];
Grav_out=[Blockmean_llf(:,1),Blockmean_llf(:,2),H_Function(Blockmean_llf(:,1),Blockmean_llf(:,2)),Blockmean_llf(:,3)+0.0419*2.67*H_Function(Blockmean_llf(:,1),Blockmean_llf(:,2)),Blockmean_Err,Blockmax_flags];

disp('Saving output')

save('Data/processedData/Terrestrial_Gravity.mat','Grav_out')

disp('Done!!')

% figure(1)
% clf
% imagesc(FreeAir_Errorm)
% block avarge bouegure to avoid alising


% figure(1)
% plot(Blockmax_flags,'.')
% title('Blockmax')
% 
% figure(2)
% plot(Blockmin_flags,'.')
% title('Blockmin')
% 
% figure(3)
% plot(Blockmean_flags,'.')
% title('Blockmean')
% 
% figure(4)
% plot(Blockmedian_flags,'.')
% title('Blockmedian')
% 
% % number of terresterial flags 
% sum(Blockmin_flags(:) == 1);%799031
% sum(Blockmax_flags(:) == 1);%795960