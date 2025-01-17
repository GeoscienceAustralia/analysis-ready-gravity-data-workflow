close all
clear 

filename='Data/GRAVITY_GRAD/Xcalibur_FVD_GDD.mat';

GravityGradient5D=importdata(filename);

% Plot 1: Gravity Gradient NSW
figure;
scatter(GravityGradient5D(:, 1), GravityGradient5D(:, 2), 1, GravityGradient5D(:, 4));
colorbar;
colormap(jet);
title(colorbar, 'mGal/m', 'FontSize', 10);
title('Gravity Gradient NSW');
xlabel('Longitude');
ylabel('Latitude');
grid on;
saveas(gcf, 'GravityGradientNSW.png'); % Save the figure
disp('Figure "GravityGradientNSW.png" saved.');

% Plot 2: Height NSW
figure;
scatter(GravityGradient5D(:, 1), GravityGradient5D(:, 2), 1, GravityGradient5D(:, 3));
colorbar;
colormap(jet);
title(colorbar, 'm', 'FontSize', 10);
title('Height NSW');
xlabel('Longitude');
ylabel('Latitude');
grid on;
saveas(gcf, 'HeightNSW.png'); % Save the figure
disp('Figure "HeightNSW.png" saved.');

% Plot 3: Error NSW
figure;
scatter(GravityGradient5D(:, 1), GravityGradient5D(:, 2), 1, GravityGradient5D(:, 5));
colorbar;
colormap(jet);
title(colorbar, 'mGal/m', 'FontSize', 10);
title('Error NSW');
xlabel('Longitude');
ylabel('Latitude');
grid on;
saveas(gcf, 'ErrorNSW.png'); % Save the figure
disp('Figure "ErrorNSW.png" saved.');



filename='Data/GRAVITY_GRAD/Otway.mat';

GravityGradient5D=importdata(filename);

% Plot 1: Gravity Gradient Otway
figure;
scatter(GravityGradient5D(:, 1), GravityGradient5D(:, 2), 1, GravityGradient5D(:, 4) * 10^(-3));
colorbar;
colormap(jet);
title(colorbar, 'mGal/m', 'FontSize', 10);
title('Gravity Gradient Otway');
xlabel('Longitude');
ylabel('Latitude');
grid on;
saveas(gcf, 'GravityGradientOtway.png'); % Save the figure
disp('Figure "GravityGradientOtway.png" saved.');

% Plot 2: Height Otway
figure;
scatter(GravityGradient5D(:, 1), GravityGradient5D(:, 2), 1, GravityGradient5D(:, 3));
colorbar;
colormap(jet);
title(colorbar, 'm', 'FontSize', 10);
title('Height Otway');
xlabel('Longitude');
ylabel('Latitude');
grid on;
saveas(gcf, 'HeightOtway.png'); % Save the figure
disp('Figure "HeightOtway.png" saved.');

% Plot 3: Error Otway
figure;
scatter(GravityGradient5D(:, 1), GravityGradient5D(:, 2), 1, GravityGradient5D(:, 5)* 10^(-3));
colorbar;
colormap(jet);
title(colorbar, 'mGal/m', 'FontSize', 10);
title('Error Otway');
xlabel('Longitude');
ylabel('Latitude');
grid on;
saveas(gcf, 'ErrorOtway.png'); % Save the figure
disp('Figure "ErrorOtway.png" saved.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import and process data
Vic_Data = importdata('Data/GRAVITY/AIRBORNE/Gippsland/GRAV.DAT');
Vic_Long = round(Vic_Data(:, 11) * 60) / 60;
Vic_Lat = round(Vic_Data(:, 10) * 60) / 60;
Vic_H = Vic_Data(:, 13);
Vic_Grav_anom = Vic_Data(:, 24);

% Plot 1: Gravity Anomaly
figure;
scatter(Vic_Long, Vic_Lat, 1, Vic_Grav_anom);
colorbar;
colormap(jet);
title(colorbar, 'mGal', 'FontSize', 10);
title('Gippsland: Gravity Anomaly');
xlabel('Longitude');
ylabel('Latitude');
grid on;
saveas(gcf, 'Gippsland_GravityAnomaly.png'); % Save the figure
disp('Figure "Gippsland_GravityAnomaly.png" saved.');

% Plot 2: Height
figure;
scatter(Vic_Long, Vic_Lat, 1, Vic_H);
colorbar;
colormap(jet);
title(colorbar, 'm', 'FontSize', 10);
title('Gippsland: Height');
xlabel('Longitude');
ylabel('Latitude');
grid on;
saveas(gcf, 'Gippsland_Height.png'); % Save the figure
disp('Figure "Gippsland_Height.png" saved.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AirborneGravitymat=importdata('Data\GRAVITY\AIRBORNE/23102024victoriaOtter/Airborne_Gravity.mat');

figure
scatter(AirborneGravitymat(:,1),AirborneGravitymat(:,2),1,AirborneGravitymat(:,4))
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('otterMat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vic_Data=importdata('Data\GRAVITY\AIRBORNE/23102024victoriaOtter/FD012_Grav.csv');
Vic_Data.data(end,:)=[];
Vic_Long=round(Vic_Data.data(:,11)*60)/60;
Vic_Lat=round(Vic_Data.data(:,10)*60)/60;
Vic_H=Vic_Data.data(:,9);
Vic_Grav_anom=Vic_Data.data(:,36);

figure
scatter(Vic_Long,Vic_Lat,1,Vic_Grav_anom)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('otterCsv')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vic_Data=importdata('Data\GRAVITY\AIRBORNE/24102024victoriaCaravan/dlv009.csv');
Vic_Data.data(end,:)=[];
Vic_Long=round(Vic_Data.data(:,35)*60)/60;
Vic_Lat=round(Vic_Data.data(:,32)*60)/60;
Vic_H=Vic_Data.data(:,39);
Vic_Grav_anom=Vic_Data.data(:,18);

figure
scatter(Vic_Long,Vic_Lat,1,Vic_Grav_anom)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('caravanCsv')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vic_Data=importdata('Data\GRAVITY\AIRBORNE/18042024victoriaNSW/dlv024_gravVictoria.csv');
Vic_Data.data(end,:)=[];
Vic_Long=round(Vic_Data.data(:,11)*60)/60;
Vic_Lat=round(Vic_Data.data(:,10)*60)/60;
Vic_H=Vic_Data.data(:,9);
Vic_Grav_anom=Vic_Data.data(:,36);

figure
scatter(Vic_Long,Vic_Lat,1,Vic_Grav_anom)
colorbar
colormap(jet)
title(colorbar,'mGal','FontSize',10);
title('aprilCsv')
