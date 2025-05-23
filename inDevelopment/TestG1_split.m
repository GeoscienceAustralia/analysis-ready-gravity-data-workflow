close all
clear all
%% Import datasets
AUDEM=importdata('AUSDEM1m.xyz');
FreeAir=importdata('FreeAir.xyz');
H=reshape(AUDEM(:,3),4861,3181)';
FA=reshape(FreeAir(:,3),4861,3181)';
Long=reshape(AUDEM(:,1),4861,3181)';
Lat=reshape(AUDEM(:,2),4861,3181)';

G1_o=FFTG1Deg(H,FA,Lat,mean(Lat(:,1)),Long,1/60,1);

ResGeoid=Stokes(Lat,Long,G1_o);

[G1_p1,G1_p2]=FFTG1Deg_split(H,FA,Lat,mean(Lat(:,1)),Long,1/60,1);

% ResGeoid_2p=Stokes(Lat,Long,G1_p1+H.*G1_p2);
% ResGeoid_3p=Stokes(Lat,Long,G1_p1)+H.*Stokes(Lat,Long,G1_p2);

figure(1)
clf
imagesc((G1_o))
colorbar
colormap(jet)
caxis([-2.5 2.5])

figure(2)
clf
imagesc((G1_p1+H.*G1_p2-G1_o))
colorbar
colormap(jet)

% 
% figure(3)
% clf
% imagesc(ResGeoid)
% colorbar
% colormap(jet)
% caxis([-0.5*10^-4 0.5*10^-4])
% 
% figure(4)
% clf
% imagesc(ResGeoid_2p)
% colorbar
% colormap(jet)
% caxis([-10^-4 10^-4])
% 
% figure(5)
% clf
% imagesc(ResGeoid_3p)
% colorbar
% colormap(jet)
% caxis([-10^-4 10^-4])


fid=fopen(['G1_0.ASC'],'w');
fprintf(fid,' %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f \n',...
    [Long(:),Lat(:),G1_o(:)]');
fclose all

fid=fopen(['G1_p1.ASC'],'w');
fprintf(fid,' %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f \n',...
    [Long(:),Lat(:),G1_p1(:)]');
fclose all

fid=fopen(['G1_p2.ASC'],'w');
fprintf(fid,' %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f %6.9f \n',...
    [Long(:),Lat(:),G1_p2(:)]');
fclose all


