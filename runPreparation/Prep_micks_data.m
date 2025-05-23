close all
clear all

MF=importdata('Lev_MF.mat');

Lev_AHD=[MF(:,1),MF(:,2),MF(:,3)-MF(:,4)];

Lev_NOC=[MF(:,1),MF(:,2),MF(:,3)-MF(:,5)];

Lev_NC=[MF(:,1),MF(:,2),MF(:,3)-MF(:,6)];

Lev_CARS=[MF(:,1),MF(:,2),MF(:,3)-MF(:,7)];

figure(1)
clf
scatter(Lev_NOC(:,1),Lev_NOC(:,2),10,Lev_NOC(:,3))
colorbar

figure(2)
clf
scatter(Lev_NOC(:,1),Lev_NOC(:,2),10,Lev_NC(:,3)-Lev_NOC(:,3))
colorbar

figure(3)
clf
scatter(Lev_NOC(:,1),Lev_NOC(:,2),10,Lev_CARS(:,3)-Lev_AHD(:,3))
colorbar
caxis([0.1 0.2])