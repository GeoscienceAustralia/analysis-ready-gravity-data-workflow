% 1.	Pilbara, Northern Western Australia (https://ecat.ga.gov.au/geonetwork/srv/api/records/33fc8e01-77d4-4a21-a8fb-452f847ac5e8)
% 2.	Eastern part of NSW (https://minview.geoscience.nsw.gov.au/#/?lon=148.5&lat=-32.5&z=7&l=)
close all
clear 
% Add the path to the function files.
addpath('functions');
% Look at NSW data.
filenameV = 'Data\airborneVectorTest\AIR_2024_Sander_Airborne_Grav_1027\AIR_2024_Sander_Airborne_Grav_1027\Located_Data\Original\Raw\Vertical_Grav_Raw.DAT';  

NSW_VData=importdata(filenameV);
NSW_VData.data(end,:)=[];
NSW_Long=round(NSW_VData.data(:,10)*60)/60;
NSW_Lat=round(NSW_VData.data(:,11)*60)/60;
NSW_H=NSW_VData.data(:,15);
NSW_VGrav=NSW_VData.data(:,39)/10;% its in mico m's per secon ^2 for some reason. i.e. needs dividing by 10 to get into mGal.

% figure
% scatter(NSW_Long, NSW_Lat, 1, NSW_H)
% colormap(jet)
% cb1 = colorbar;                       
% title(cb1,'m','FontSize',10)
% title('Flight height')
% saveas(gcf,'outputs/plots/NSWheight.png');
% 
% figure
% scatter(NSW_Long, NSW_Lat, 1, NSW_VGrav)
% colormap(jet)
% cb1 = colorbar;                      
% title(cb1,'mGal','FontSize',10)
% title('Vertical')
% saveas(gcf,'outputs/plots/NSWvertical.png');


filenameH = 'Data\airborneVectorTest\AIR_2024_Sander_Airborne_Grav_1027\AIR_2024_Sander_Airborne_Grav_1027\Located_Data\Original\Raw\Horizontal_Grav_Raw.DAT';

NSW_HData=importdata(filenameH);
NSW_HData.data(end,:)=[];
NSW_EGrav=NSW_HData.data(:,32)/10;% its in mico m's per secon ^2 for some reason. i.e. needs dividing by 10 to get into mGal.
NSW_NGrav=NSW_HData.data(:,43)/10;% its in mico m's per secon ^2 for some reason. i.e. needs dividing by 10 to get into mGal.
NSW_EGrav(NSW_EGrav < -1000) = NaN;
NSW_NGrav(NSW_NGrav < -1000) = NaN;

% figure
% scatter(NSW_Long, NSW_Lat, 1, NSW_EGrav)
% colormap(jet)
% cb1 = colorbar;                      
% title(cb1,'mGal','FontSize',10)
% title('East')
% saveas(gcf,'outputs/plots/NSWeast.png');
% 
% figure
% scatter(NSW_Long, NSW_Lat, 1, NSW_NGrav)
% colormap(jet)
% cb1 = colorbar;                       
% title(cb1,'mGal','FontSize',10)
% title('North')
% saveas(gcf,'outputs/plots/NSWnorth.png');

disp('Flight height')
displayStats(NSW_H)
disp('East')
displayStats(NSW_EGrav)
disp('North')
displayStats(NSW_NGrav)
disp('Vertical')
displayStats(NSW_VGrav)

% Combine data sets. East and North has been added.
ABGrav=[NSW_Long, NSW_Lat, NSW_H, NSW_VGrav, NSW_VGrav*0+3, NSW_EGrav, NSW_NGrav];

% % Block mean as other datasets.
AB_Grav_Long=round(ABGrav(:,1)*60)/60;   
AB_Grav_Lat=round(ABGrav(:,2)*60)/60;
[vals, uid, idx]=unique([AB_Grav_Long,AB_Grav_Lat],'rows');
AB_Grav_BM(:,3) = accumarray(idx,ABGrav(:,3),[],@mean);
AB_Grav_BM(:,4) = accumarray(idx,ABGrav(:,4),[],@mean);
AB_Grav_BM(:,5) = accumarray(idx,ABGrav(:,5),[],@mean);
AB_Grav_BM(:,7) = accumarray(idx,ABGrav(:,6),[],@mean);% East  added.
AB_Grav_BM(:,8) = accumarray(idx,ABGrav(:,7),[],@mean);% North added.
AB_Grav_BM(:,1) = AB_Grav_Long(uid);
AB_Grav_BM(:,2) = AB_Grav_Lat(uid);

AB_Grav_BM(:,6)=3;% Flag for airborne gravity
save('Data\processedData\Airborne3DNSW.mat','AB_Grav_BM')






% 
% subplot(1,4,2)
% hold on
% scatter(NSW_Long, NSW_Lat, 1, NSW_Data.data(:,76))
% colormap(jet)
% cb1 = colorbar;                       % capture handle so we can title it
% title(cb1,'mGal','FontSize',10)
% title('East')
% 
% subplot(1,4,3)
% hold on
% scatter(NSW_Long, NSW_Lat, 1, NSW_Data.data(:,77))
% colormap(jet)
% cb1 = colorbar;                       % capture handle so we can title it
% title(cb1,'mGal','FontSize',10)
% title('North')
% 
% subplot(1,4,4)
% hold on
% scatter(NSW_Long, NSW_Lat, 1, NSW_Grav_anom)
% colormap(jet)
% cb1 = colorbar;                       % capture handle so we can title it
% title(cb1,'mGal','FontSize',10)
% title('Vertical')




% Check by comparing to EGM2008
% GGM=importdata('Data/GGM/EGM2008_For_Gridded_Int.mat');
% GGM_Gi=griddedInterpolant(GGM.x,GGM.y,GGM.z,GGM.g);
% 
% GGM_Gi_interpolatedNSW=GGM_Gi(NSW_Long,-NSW_Lat,NSW_H);



% figure
% 
% % ---------- Sub‑plot 1 ----------
% subplot(2,1,1)
% hold on
% scatter(NSW_Long, NSW_Lat, 1, NSW_Grav_anom)
% colormap(jet)
% cb1 = colorbar;                       % capture handle so we can title it
% title(cb1,'mGal','FontSize',10)
% title('Vertical')
% 
% % ---------- Sub‑plot 2 ----------
% subplot(2,1,2)
% hold on
% diffVals = NSW_Grav_anom - GGM_Gi_interpolatedNSW;
% scatter(NSW_Long, NSW_Lat, 1, diffVals)
% colormap(jet)
% cb2 = colorbar;
% title(cb2,'mGal','FontSize',10)
% title('Vertical - EGM2008')
% 
% % ---------- Compute & display mean ----------
% meanDiff = mean(diffVals);
% % ---------- Add mean as figure text ----------
% % Place it near the bottom of the figure; adjust position to taste
% annotation('textbox',[0.15 0.02 0.7 0.05], ...
%            'String',sprintf('Mean difference (GGM – Vertical): %.4f mGal',meanDiff), ...
%            'EdgeColor','none', ...
%            'HorizontalAlignment','center', ...
%            'FontSize',10, ...
%            'FontWeight','bold');

% ---------- Save ----------
%saveas(gcf,'outputs/plots/EGM2008Victoria.png');
