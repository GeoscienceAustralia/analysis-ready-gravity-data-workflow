close all
clear 

% Look at Victoria Data first.
Vic_Data=importdata('Data\GRAVITY\AIRBORNE/23102024victoriaOtter/FD012_Grav.csv');
Vic_Data.data(end,:)=[];
Vic_LongOtter=round(Vic_Data.data(:,11)*60)/60;
Vic_LatOtter=round(Vic_Data.data(:,10)*60)/60;
Vic_HOtter=Vic_Data.data(:,9);
Vic_Grav_anomOtter=Vic_Data.data(:,36);

% this file has different columns
Vic_Data=importdata('Data\GRAVITY\AIRBORNE/24102024victoriaCaravan/dlv009.csv');
Vic_Data.data(end,:)=[];
Vic_LongCaravan=round(Vic_Data.data(:,35)*60)/60;
Vic_LatCaravan=round(Vic_Data.data(:,32)*60)/60;
Vic_HCaravan=Vic_Data.data(:,39);
Vic_Grav_anomCaravan=Vic_Data.data(:,18);

figure
scatter(Vic_LongOtter,Vic_LatOtter,1,Vic_Grav_anomOtter)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('Otter')

figure
scatter(Vic_LongCaravan,Vic_LatCaravan,1,Vic_Grav_anomCaravan)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('Caravan')


figure
scatter(Vic_LongOtter,Vic_LatOtter,1,Vic_Grav_anomOtter)
hold on
scatter(Vic_LongCaravan,Vic_LatCaravan,1,Vic_Grav_anomCaravan)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('OtterCaravan')

% Check by comparing to EGM2008
GGM=importdata('Data/GGM/EGM2008_For_Gridded_Int.mat');
GGM_Gi=griddedInterpolant(GGM.x,GGM.y,GGM.z,GGM.g);
GGM_Gi_interpolatedCaravan=GGM_Gi(Vic_LongCaravan,-Vic_LatCaravan,Vic_HCaravan);
GGM_Gi_interpolatedOtter=GGM_Gi(Vic_LongOtter,-Vic_LatOtter,Vic_HOtter);

figure
subplot(2, 1, 1) % Create a subplot with 1 row and 2 columns, and set the first subplot as active
hold on
scatter(Vic_LongOtter,Vic_LatOtter,1,Vic_Grav_anomOtter)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('Otter')

subplot(2, 1, 2) % Create a subplot with 1 row and 2 columns, and set the first subplot as active
hold on
scatter(Vic_LongOtter,Vic_LatOtter,1,Vic_Grav_anomOtter-GGM_Gi_interpolatedOtter)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('Otter-EGM2008')

disp('Mean of Difference GGM and Otter ')
mean(Vic_Grav_anomOtter-GGM_Gi_interpolatedOtter)

figure
subplot(2, 1, 1) % Create a subplot with 1 row and 2 columns, and set the first subplot as active
hold on
scatter(Vic_LongCaravan,Vic_LatCaravan,1,Vic_Grav_anomCaravan)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('Caravan')

subplot(2, 1, 2) % Create a subplot with 1 row and 2 columns, and set the first subplot as active
hold on
scatter(Vic_LongCaravan,Vic_LatCaravan,1,Vic_Grav_anomCaravan-GGM_Gi_interpolatedCaravan)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('Caravan-EGM2008')

disp('Mean of Difference GGM and Caravan ')
mean(Vic_Grav_anomCaravan-GGM_Gi_interpolatedCaravan)

% Look at NSW second.

NSW_Data=importdata('Data\GRAVITY\AIRBORNE\18042024victoriaNSW/gravNSW.csv');
NSW_Data.data(end,:)=[];
NSW_Long=round(NSW_Data.data(:,65)*60)/60;
NSW_Lat=round(NSW_Data.data(:,62)*60)/60;
NSW_H=NSW_Data.data(:,66);
NSW_Grav_anom=NSW_Data.data(:,36)/10;% its in mico m's per secon ^2 for some reason. i.e. needs dividing by 10 to get into mGal.

GGM_Gi_interpolatedNSW=GGM_Gi(NSW_Long,-NSW_Lat,NSW_H);

figure
subplot(2, 1, 1) % Create a subplot with 1 row and 2 columns, and set the first subplot as active
hold on
scatter(NSW_Long,NSW_Lat,1,NSW_Grav_anom)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('NSW')

subplot(2, 1, 2) % Create a subplot with 1 row and 2 columns, and set the first subplot as active
hold on
scatter(NSW_Long,NSW_Lat,1,NSW_Grav_anom-GGM_Gi_interpolatedNSW)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('NSW-EGM2008')

disp('Mean of Difference GGM and NSW ')
mean(NSW_Grav_anom-GGM_Gi_interpolatedNSW)

% Combine data sets.
ABGrav=[Vic_LongOtter,Vic_LatOtter,Vic_HOtter,Vic_Grav_anomOtter,Vic_Grav_anomOtter*0+3;...
        Vic_LongCaravan,Vic_LatCaravan,Vic_HCaravan,Vic_Grav_anomCaravan,Vic_Grav_anomCaravan*0+3;...
        NSW_Long,NSW_Lat,NSW_H,NSW_Grav_anom,NSW_Grav_anom*0+3];

% % Block mean as other datasets.
AB_Grav_Long=round(ABGrav(:,1)*60)/60;   
AB_Grav_Lat=round(ABGrav(:,2)*60)/60;
[vals, uid, idx]=unique([AB_Grav_Long,AB_Grav_Lat],'rows');
AB_Grav_BM(:,3) = accumarray(idx,ABGrav(:,3),[],@mean);
AB_Grav_BM(:,4) = accumarray(idx,ABGrav(:,4),[],@mean);
AB_Grav_BM(:,5) = accumarray(idx,ABGrav(:,5),[],@mean);
AB_Grav_BM(:,1) = AB_Grav_Long(uid);
AB_Grav_BM(:,2) = AB_Grav_Lat(uid);
% 
save('Data\processedData\AirborneGravityCaravanOtterNSW.mat','AB_Grav_BM')
% disp('Done')


