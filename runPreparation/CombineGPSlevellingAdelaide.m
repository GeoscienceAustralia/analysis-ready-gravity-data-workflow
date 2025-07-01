% combine gravity gradient data
Otway=importdata('Data/GRAVITY_GRAD/OtwayMgalm.mat');
Xcalibur=importdata('Data/GRAVITY_GRAD/Xcalibur_FVD_GDD.mat');

% Combine data sets.
all2=[Xcalibur;Otway];

save('Data/GRAVITY_GRAD/XcaliburOtway.mat','all2')


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

save('Data/GPS_LEVELLING/Lev_Adelaide.mat','all3GPSlevelling')


