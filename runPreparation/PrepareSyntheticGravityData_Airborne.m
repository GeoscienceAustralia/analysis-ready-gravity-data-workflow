close all
clear 

% Look at Data first.
Data=importdata('Data\GRAVITY\AIRBORNE/18122024PerthSynthetic/Gravity_500mDrape_plusNoise_3mGal_70s_LP_1kLines.csv');
%Data1=importdata('Data\GRAVITY\AIRBORNE/18122024PerthSynthetic/InputGravity_500mDrape_2kLines.csv');
%Data3=importdata('Data\GRAVITY\AIRBORNE/13122024PerthSynthetic/SyntheticGravity_Perth.csv');

% figure
% scatter(Data3.data(:,10), Data3.data(:,9), 1, Data3.data(:,12))
% colorbar
% colormap(jet)
% title(colorbar, 'mGal', 'FontSize', 10);
% title('Height ')
% 
% 
% figure
% scatter(Data3.data(:,10), Data3.data(:,9), 1, Data3.data(:,1))
% colorbar
% colormap(jet)
% title(colorbar, 'mGal', 'FontSize', 10);
% title('Height ')

Data.data(end,:)=[];
Long=round(Data.data(:,3)*60)/60;
Lat=round(Data.data(:,2)*60)/60;
Height=Data.data(:,6);
Grav_anom=Data.data(:,1);

% Check by comparing to EGM2008
GGM=importdata('Data/GGM/EGM2008_For_Gridded_Int.mat');
GGM_Gi=griddedInterpolant(GGM.x,GGM.y,GGM.z,GGM.g);
GGM_Gi_interpolated=GGM_Gi(Long,-Lat,Height);
% 

figure
scatter(Long, Lat, 1, Height)
colorbar
colormap(jet)
title(colorbar, 'mGal', 'FontSize', 10);
title('Height ')

figure
subplot(1, 2, 1) % Create a subplot with 1 row and 2 columns, and set the first subplot as active
hold on
scatter(Long, Lat, 1, Grav_anom)
colorbar
colormap(jet)
title(colorbar, 'mGal', 'FontSize', 10);
title('airborneGravity')

subplot(1, 2, 2) % Set the second subplot as active
hold on
scatter(Long, Lat, 1, Grav_anom - GGM_Gi_interpolated)
colorbar
colormap(jet)
title(colorbar, 'mGal', 'FontSize', 10);
title('airborneGravity-EGM2008')

disp('Mean of Difference GGM and Airborne ')
mean(Grav_anom - GGM_Gi_interpolated)

% Combine data sets.
ABGrav=[Long,Lat,Height,Grav_anom,Grav_anom*0+3];
% Block mean as other datasets.
AB_Grav_Long=round(ABGrav(:,1)*60)/60;   
AB_Grav_Lat=round(ABGrav(:,2)*60)/60;
[vals, uid, idx]=unique([AB_Grav_Long,AB_Grav_Lat],'rows');
AB_Grav_BM(:,3) = accumarray(idx,ABGrav(:,3),[],@mean);
AB_Grav_BM(:,4) = accumarray(idx,ABGrav(:,4),[],@mean);
AB_Grav_BM(:,5) = accumarray(idx,ABGrav(:,5),[],@mean);
AB_Grav_BM(:,1) = AB_Grav_Long(uid);
AB_Grav_BM(:,2) = AB_Grav_Lat(uid);

save('Data\processedData\AirborneGravityPerthSynthetic70sLP1kLines.mat','AB_Grav_BM')
disp('Done')





























