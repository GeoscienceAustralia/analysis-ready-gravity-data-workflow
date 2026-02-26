% Specify file name and sheet (optional)
% filename = 'Data/GPS_LEVELLING/taje_a_1412353_sm4972.xls';
% sheetName ='AHD'; % 'CARS2009'; %'AHD'; 
% 
% % Read the data
% data = readtable(filename, 'Sheet', sheetName);
% 
% % Keep only rows up to 7000 (or the max number of rows if fewer exist)
% maxRow = min(7319, height(data));
% data = data(1:maxRow, :);
% 
% % Extract three specific columns by name
% %col0 = data.IDNO;
% col1 = data.LONDD;
% col2 = data.LATDD;
% col3 = data.Zeta;
% 
% % Combine into a single array if needed
% selectedData = [col1 col2 col3];
% 
% save('Data/GPS_LEVELLING/AHDzeta7319.mat','selectedData')
%save('Data/GPS_LEVELLING/CARS2009zeta7301.mat','selectedData')



Lev131=importdata('Data/GPS_LEVELLING/Lev131pointsNEXY.mat');
Lev131(Lev131(:,3) == 1000, :) = [];
Lev2=importdata('Data/GPS_LEVELLING/Lev2_97pointsAdelaideHills.mat');
Lev3=importdata('Data/GPS_LEVELLING/Lev3_92pointsSouthNEXY.mat');


% common variables for plotting
axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;
% 
% 
scatter(Lev131(:,1),Lev131(:,2),15,Lev131(:,3)*0+1,'filled')
hold on 
scatter(Lev2(:,1),Lev2(:,2),15,Lev2(:,3)*0+2,'filled')
hold on 
scatter(Lev3(:,1),Lev3(:,2),15,Lev3(:,3)*0+3,'filled')
customizeMap('',' ',Coastline,axisLimits)

% combine GPS-levelling data for greater Adelaide  
figure
scatter(Lev2(:,1),Lev2(:,2),10,Lev2(:,3))
colorbar
colormap(jet)
title(colorbar,'m','FontSize',10);
title('Lev2')

figure
scatter(Lev3(:,1),Lev3(:,2),10,Lev3(:,3))
colorbar
colormap(jet)
title(colorbar,'m','FontSize',10);
title('Lev3')

Lev(Lev(:,3) == 1000, :) = [];

% Combine data sets.
all3GPSlevelling=[Lev(:,1),Lev(:,2),Lev(:,3);...
                  Lev2(:,1),Lev2(:,2),Lev2(:,3);...
                  Lev3(:,1),Lev3(:,2),Lev3(:,3)];

figure
scatter(all3GPSlevelling(:,1),all3GPSlevelling(:,2),10,all3GPSlevelling(:,3))
colorbar
colormap(jet)
title(colorbar,'m','FontSize',10);
% 
% save('Data/GPS_LEVELLING/Lev_Adelaide.mat','all3GPSlevelling')








