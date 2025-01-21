close all
clear 

eotvos2mgalm = 1e-4;

filename = 'Data/GRAVITY_GRAD/Otway.mat';

GravityGradient5DOtway = importdata(filename);

% Multiply columns 4 and 5 by eotvos2mgalm
GravityGradient5DOtway(:, 4) = GravityGradient5DOtway(:, 4) * eotvos2mgalm;
GravityGradient5DOtway(:, 5) = GravityGradient5DOtway(:, 5) * eotvos2mgalm;

% Save the modified data to a .mat file
save('OtwayMgalm.mat', 'GravityGradient5DOtway');

% filename='Data/GRAVITY_GRAD/Xcalibur_FVD_GDD.mat';
% 
% GravityGradient5DNSW=importdata(filename);

% Plot 1: Gravity Gradient NSW
figure;
scatter(GravityGradient5DNSW(:, 1), GravityGradient5DNSW(:, 2), 1, GravityGradient5DNSW(:, 4));
hold on
scatter(GravityGradient5DOtway(:, 1), GravityGradient5DOtway(:, 2), 1, GravityGradient5DOtway(:, 4)* eotvos2mgalm);
colorbar;
colormap(jet);
title(colorbar, 'mGal/m', 'FontSize', 10);
title('Gravity Gradient');
xlabel('Longitude');
ylabel('Latitude');
grid on;
saveas(gcf, 'GravityGradient.png'); % Save the figure

% 
% % Plot 2: Height NSW
% figure;
% scatter(GravityGradient5DNSW(:, 1), GravityGradient5DNSW(:, 2), 1, GravityGradient5DNSW(:, 3));
% hold on
% scatter(GravityGradient5DOtway(:, 1), GravityGradient5DOtway(:, 2), 1, GravityGradient5DOtway(:, 3));
% colorbar;
% colormap(jet);
% title(colorbar, 'm', 'FontSize', 10);
% title('Height');
% xlabel('Longitude');
% ylabel('Latitude');
% grid on;
% saveas(gcf, 'Height.png'); % Save the figure

% 
% % Plot 3: Error NSW
% figure;
% scatter(GravityGradient5DNSW(:, 1), GravityGradient5DNSW(:, 2), 1, GravityGradient5DNSW(:, 5));
% hold on
% scatter(GravityGradient5DOtway(:, 1), GravityGradient5DOtway(:, 2), 1, GravityGradient5DOtway(:, 5)* eotvos2mgalm);
% colorbar;
% colormap(jet);
% title(colorbar, 'mGal/m', 'FontSize', 10);
% title('Error');
% xlabel('Longitude');
% ylabel('Latitude');
% grid on;
% saveas(gcf, 'Error.png'); % Save the figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AirborneGippsland = importdata('Data/GRAVITY/AIRBORNE/Gippsland/GRAV.DAT');
AirborneGravitymat=importdata('Data\GRAVITY\AIRBORNE/23102024victoriaOtter/Airborne_Gravity.mat');

% Plot 1: Gravity Anomaly
figure;
scatter(AirborneGippsland(:,10),AirborneGippsland(:,11),1,AirborneGippsland(:,36));
hold on
scatter(AirborneGravitymat(:,1),AirborneGravitymat(:,2),1,AirborneGravitymat(:,4))
colorbar;
colormap(jet);
title(colorbar, 'mGal', 'FontSize', 10);
title('Gravity Anomaly');
xlabel('Longitude');
ylabel('Latitude');
grid on;
%saveas(gcf, 'GippslandGravityAnomaly.png'); % Save the figure

% Plot 2: Height
figure;
scatter(AirborneGippsland(:,10),AirborneGippsland(:,11),1,AirborneGippsland(:,13));
hold on
scatter(AirborneGravitymat(:,1),AirborneGravitymat(:,2),1,AirborneGravitymat(:,3))
colorbar;
colormap(jet);
title(colorbar, 'm', 'FontSize', 10);
title('Height');
xlabel('Longitude');
ylabel('Latitude');
grid on;

% Plot 2: Height
figure;
scatter(AirborneGippsland(:,10),AirborneGippsland(:,11),1,'b');
% hold on
%scatter(AirborneGravitymat(:,1),AirborneGravitymat(:,2),1,'r')
colorbar;
colormap(jet);
title(colorbar, 'm', 'FontSize', 10);
title('Height');
xlabel('Longitude');
ylabel('Latitude');
grid on;

%saveas(gcf, 'GippslandHeight.png'); % Save the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vic_Data=importdata('Data\GRAVITY\AIRBORNE/23102024victoriaOtter/FD012_Grav.csv');
Vic_Data.data(end,:)=[];
Vic_Long=round(Vic_Data.data(:,11)*60)/60;
Vic_Lat=round(Vic_Data.data(:,10)*60)/60;
Vic_H=Vic_Data.data(:,9);
Vic_Grav_anom=Vic_Data.data(:,36);

figure
scatter(Vic_Long,Vic_Lat,1,Vic_Grav_anom)
hold on
scatter(AirborneGravitymat(:,1),AirborneGravitymat(:,2),1,AirborneGravitymat(:,4))
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('otterCsv')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vic_Data=importdata('Data\GRAVITY\AIRBORNE/24102024victoriaCaravan/dlv009.csv');
% Vic_Data.data(end,:)=[];
% Vic_Long=round(Vic_Data.data(:,35)*60)/60;
% Vic_Lat=round(Vic_Data.data(:,32)*60)/60;
% Vic_H=Vic_Data.data(:,39);
% Vic_Grav_anom=Vic_Data.data(:,18);
% 
% figure
% scatter(Vic_Long,Vic_Lat,1,Vic_Grav_anom)
% colorbar
% colormap(jet)
% title(colorbar,'mGal','FontSize',10);
% title('caravanCsv')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vic_Data=importdata('Data\GRAVITY\AIRBORNE/18042024victoriaNSW/dlv024_gravVictoria.csv');
% Vic_Data.data(end,:)=[];
% Vic_Long=round(Vic_Data.data(:,11)*60)/60;
% Vic_Lat=round(Vic_Data.data(:,10)*60)/60;
% Vic_H=Vic_Data.data(:,9);
% Vic_Grav_anom=Vic_Data.data(:,36);
% 
% figure
% scatter(Vic_Long,Vic_Lat,1,Vic_Grav_anom)
% colorbar
% colormap(jet)
% title(colorbar,'mGal','FontSize',10);
% title('aprilCsv')



