close all
clear 
%% Import the various gravity datasets
Grav1=importdata('Data\processedData\AirborneAllJuly14.mat');
%Grav1=importdata('Data\processedData\AirborneAll.mat');
%Grav1=importdata('Data\processedData\AirborneGravityGippslandCaravanOtterNSW.mat');
%Grav1=importdata('Data\processedData\AirborneGravityCaravanOtterNSW.mat');
%Grav1=importdata('Data\processedData\AirborneGravityPerthSynthetic70sLP1kLines.mat');
%Grav1=importdata('Data\processedData\AirborneGravityVicNSW.mat');
%Grav2=importdata('AIRBORNE/WA/WA_Airborne.mat');
%Grav3=importdata('AIRBORNE/NSW/Airborne_Gravity_SGL.mat');% look at more scripts here, there is bias here.
% relative to EGM2008, calculating bias at intersection.
% Grav3=importdata('AIRBORNE/NSW/Xcalibur_Gravity.mat');
Grav4=importdata('Data\processedData\Terrestrial_Gravity.mat');
Grav1(:,6)=3;% Flag for airborne gravity
%% Create output
% Grav_out=[Longitude,Latitude,Ortho-H,gravity anomaly,uncertainty,flag];
% Out=[Grav1;Grav2;Grav3;Grav4];GGM_Zetai
%Out=[Grav1;Grav2;Grav3;Grav4];
Out=[Grav1;Grav4];
Out(isnan(Out(:,5)),:)=[];
Out(isnan(Out(:,4)),:)=[];
Out(isnan(Out(:,3)),:)=[];
%save('Data\processedData\GravityAllVicNSW.mat','Out')
%save('Data\processedData\GravityAllPerthSynthetic70sLP1kLines.mat','Out')
%save('Data\processedData\GravityAllCaravanOtterNSW.mat','Out')
%save('Data\processedData\GravityAllGippslandCaravanOtterNSW.mat','Out')
save('Data\processedData\GravityAllTerrestrialAirborneJuly14.mat','Out')
disp('Done')