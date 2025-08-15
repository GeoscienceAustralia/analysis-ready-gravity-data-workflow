% Specify file name and sheet (optional)
filename = 'Data/GPS_LEVELLING/taje_a_1412353_sm4972.xls';
sheetName = 'CARS2009';%'AHD'; 

% Read the data
data = readtable(filename, 'Sheet', sheetName);

% Extract three specific columns by name
% Replace 'Column1', 'Column2', 'Column3' with your actual column headers
col1 = data.LONDD;
col2 = data.LATDD;
col3 = data.Zeta;

% Combine into a single array if needed
selectedData = [col1 col2 col3];

%save('Data/GPS_LEVELLING/AHDzeta7506.mat','selectedData')
save('Data/GPS_LEVELLING/CARS2009zeta7506.mat','selectedData')













% % combine GPS-levelling data for greater Adelaide  
% figure
% scatter(Lev2(:,1),Lev2(:,2),10,Lev2(:,3))
% colorbar
% colormap(jet)
% title(colorbar,'m','FontSize',10);
% title('Lev2')
% 
% figure
% scatter(Lev3(:,1),Lev3(:,2),10,Lev3(:,3))
% colorbar
% colormap(jet)
% title(colorbar,'m','FontSize',10);
% title('Lev3')
% 
% Lev(Lev(:,3) == 1000, :) = [];
% 
% % Combine data sets.
% all3GPSlevelling=[Lev(:,1),Lev(:,2),Lev(:,3);...
%                   Lev2(:,1),Lev2(:,2),Lev2(:,3);...
%                   Lev3(:,1),Lev3(:,2),Lev3(:,3)];
% 
% figure
% scatter(all3GPSlevelling(:,1),all3GPSlevelling(:,2),10,all3GPSlevelling(:,3))
% colorbar
% colormap(jet)
% title(colorbar,'m','FontSize',10);
% 
% save('Data/GPS_LEVELLING/Lev_Adelaide.mat','all3GPSlevelling')








