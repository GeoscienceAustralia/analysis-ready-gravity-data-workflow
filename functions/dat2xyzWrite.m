
% Path to AUSGeoid dat file
AUSGeoidFile ='Data/EXISTING_GEOID_MODELS/AUSGeoid2020_20180201_win.dat'; %from GA website 
% https://www.ga.gov.au/scientific-topics/positioning-navigation/positioning-australia/geodesy/ahdgm/ausgeoid2020

xyzAUSGeoid = readDatFile(AUSGeoidFile);

% Path to AGQG dat file
AGQGFile =    'Data/EXISTING_GEOID_MODELS/AGQG2017_20170907.dat';

xyzAGQG = readDatFile(AGQGFile);

%The zero degree term offset between AGQG_20201120 and AGQG_2017 is -0.41 m. 
% Heights above AGQG_2017 will be 0.41 m smaller than heights above AGQG_20201120.

xyz = [xyzAUSGeoid(:,1) xyzAUSGeoid(:,2) xyzAUSGeoid(:,3)-xyzAGQG(:,3)+0.41];

% Plot
figure;
scatter(xyz(:,1), xyz(:,2), 1, xyz(:,3), 'filled');
colorbar;
colormap(parula);
xlabel('Longitude');
ylabel('Latitude');
title('Geoid in meters');

 
%Write output file
xyzFile ='outputs/diffAUSGeoid2020AGQG2017.xyz';
writematrix(xyz, xyzFile, 'Delimiter', ' ', 'FileType', 'text');

fprintf('✅ XYZ file written: %s\n', xyzFile);
