close all
clear 

Adelaide_Data=importdata('Data/GRAVITY/AIRBORNE/adelaide2025-06-05/GRAV.DAT');



AB_Grav_BMtest=importdata('Data\GRAVITY\AIRBORNE/23102024victoriaOtter/Airborne_Gravity.mat');
figure
scatter(AB_Grav_BMtest(:,1),AB_Grav_BMtest(:,2),1,AB_Grav_BMtest(:,4))
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);

% Look at Victoria Data first.
%Vic_Data=importdata('Data\GRAVITY\AIRBORNE/18042024victoriaNSW/dlv024_gravVictoria.csv');
Vic_Data=importdata('Data\GRAVITY\AIRBORNE/23102024victoriaOtter/FD012_Grav.csv');
Vic_Data.data(end,:)=[];
Vic_Long=round(Vic_Data.data(:,11)*60)/60;
Vic_Lat=round(Vic_Data.data(:,10)*60)/60;
Vic_H=Vic_Data.data(:,9);
Vic_Grav_anom=Vic_Data.data(:,36);

% this file has different columns
% Vic_Data=importdata('Data\GRAVITY\AIRBORNE/24102024victoriaCaravan/dlv009.csv');
% Vic_Data.data(end,:)=[];
% Vic_Long=round(Vic_Data.data(:,35)*60)/60;
% Vic_Lat=round(Vic_Data.data(:,32)*60)/60;
% Vic_H=Vic_Data.data(:,39);
% Vic_Grav_anom=Vic_Data.data(:,18);

%[Vic_Latm,Vic_Longm]=meshgrid(min(Vic_Lat):1/60:max(Vic_Lat),min(Vic_Long):1/60:max(Vic_Long));

% Check by comparing to EGM2008
% GGM=importdata('Data/GGM/EGM2008_For_Gridded_Int.mat');
% GGM_Gi=griddedInterpolant(GGM.x,GGM.y,GGM.z,GGM.g);
% GGM_Gi_interpolated=GGM_Gi(Vic_Long,-Vic_Lat,Vic_H);

%figure
%subplot(2, 1, 1) % Create a subplot with 1 row and 2 columns, and set the first subplot as active
hold on
scatter(Vic_Long,Vic_Lat,1,Vic_Grav_anom)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('airborneGravity')

% subplot(2, 1, 2) % Create a subplot with 1 row and 2 columns, and set the first subplot as active
% hold on
% scatter(Vic_Long,Vic_Lat,1,Vic_Grav_anom-GGM_Gi_interpolated)
% colorbar
% colormap(jet)
% title(colorbar,'mGal','FontSize',10);
% title('airborneGravity-EGM2008')

% disp('Mean of Difference GGM and Airborne ')
% mean(Vic_Grav_anom-GGM_Gi_interpolated)
% 
% ABGrav=[Vic_Long,Vic_Lat,Vic_H,Vic_Grav_anom,Vic_Grav_anom*0+3];
% % Block mean as other datasets.
% AB_Grav_Long=round(ABGrav(:,1)*60)/60;   
% AB_Grav_Lat=round(ABGrav(:,2)*60)/60;
% [vals, uid, idx]=unique([AB_Grav_Long,AB_Grav_Lat],'rows');
% AB_Grav_BM(:,3) = accumarray(idx,ABGrav(:,3),[],@mean);
% AB_Grav_BM(:,4) = accumarray(idx,ABGrav(:,4),[],@mean);
% AB_Grav_BM(:,5) = accumarray(idx,ABGrav(:,5),[],@mean);
% AB_Grav_BM(:,1) = AB_Grav_Long(uid);
% AB_Grav_BM(:,2) = AB_Grav_Lat(uid);
% 
% save('Data\processedData\AirborneGravityVictoriaOtter.mat','AB_Grav_BM')
% disp('Done')



%% Look at NSW second.

% NSW_Data=importdata('Data\GRAVITY\AIRBORNE\18042024/grav.csv');
% NSW_Data.data(end,:)=[];
% NSW_Long=round(NSW_Data.data(:,65)*60)/60;
% NSW_Lat=round(NSW_Data.data(:,62)*60)/60;
% NSW_H=NSW_Data.data(:,66);
% NSW_Grav_anom=NSW_Data.data(:,36)/10;% its in mico m's per secon ^2 for some reason. i.e. needs dividing by 10 to get into mGal.
% 
% [NSW_Latm,NSW_Longm]=meshgrid(min(NSW_Lat):1/60:max(NSW_Lat),min(NSW_Long):1/60:max(NSW_Long));
% 
% figure(3)
% clf
% hold on
% scatter(NSW_Long,NSW_Lat,1,NSW_Grav_anom)
% colorbar
% caxis([-100 100])
% colormap(jet)
% 
% figure(4)
% clf
% hold on
% scatter(NSW_Long,NSW_Lat,1,NSW_Grav_anom-GGM_Gi(NSW_Long,-NSW_Lat,NSW_H))
% colorbar
% caxis([-100 100])
% colormap(jet)
% title('Diferences look small and zero mean? good to go!')

%% Combine data sets.
% ABGrav=[Vic_Long,Vic_Lat,Vic_H,Vic_Grav_anom,Vic_Grav_anom*0+3;...
%         NSW_Long,NSW_Lat,NSW_H,NSW_Grav_anom,NSW_Grav_anom*0+3];
% % Block mean as other datasets.
% AB_Grav_Long=round(ABGrav(:,1)*60)/60;   
% AB_Grav_Lat=round(ABGrav(:,2)*60)/60;
% [vals, uid, idx]=unique([AB_Grav_Long,AB_Grav_Lat],'rows');
% AB_Grav_BM(:,3) = accumarray(idx,ABGrav(:,3),[],@mean);
% AB_Grav_BM(:,4) = accumarray(idx,ABGrav(:,4),[],@mean);
% AB_Grav_BM(:,5) = accumarray(idx,ABGrav(:,5),[],@mean);
% AB_Grav_BM(:,1) = AB_Grav_Long(uid);
% AB_Grav_BM(:,2) = AB_Grav_Lat(uid);
% 
% save('Data\processedData\AirborneGravityVicNSW.mat','AB_Grav_BM')
% disp('Done')


