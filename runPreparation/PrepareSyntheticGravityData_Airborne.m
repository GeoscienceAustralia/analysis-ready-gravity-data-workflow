close all
clear 


% Look at Data first.
Data=importdata('Data\GRAVITY\AIRBORNE/13122024PerthSynthetic/SyntheticGravity_Perth.csv');
Data.data(end,:)=[];
Long=round(Data.data(:,10)*60)/60;
Lat=round(Data.data(:,9)*60)/60;
Height=Data.data(:,12);
Uncertain=Data.data(:,14);

% Check by comparing to EGM2008
GGM=importdata('Data/GGM/EGM2008_For_Gridded_Int.mat');
GGM_Gi=griddedInterpolant(GGM.x,GGM.y,GGM.z,GGM.g);
GGM_Gi_interpolated=GGM_Gi(Long,-Lat,Height);
% 

figure
    subplot(1, 2, 1) % Create a subplot with 1 row and 2 columns, and set the first subplot as active
    hold on
    scatter(Long, Lat, 1, Height)
    colorbar
    colormap(jet)
    title(colorbar, 'mGal', 'FontSize', 10);
    title('Height ')

    subplot(1, 2, 2) % Set the second subplot as active
    hold on
    scatter(Long, Lat, 1, Uncertain)
    colorbar
    colormap(jet)
    title(colorbar, 'mGal', 'FontSize', 10);
    title('Uncertainty ')

    % Save the figure as a PNG file
    saveas(gcf, 'UncertainyHeight.png')

for column_num = 3:8
    Grav_anom = Data.data(:, column_num);

    figure
    subplot(1, 2, 1) % Create a subplot with 1 row and 2 columns, and set the first subplot as active
    hold on
    scatter(Long, Lat, 1, Grav_anom)
    colorbar
    colormap(jet)
    title(colorbar, 'mGal', 'FontSize', 10);
    title(['column ', num2str(column_num)])

    subplot(1, 2, 2) % Set the second subplot as active
    hold on
    scatter(Long, Lat, 1, Grav_anom - GGM_Gi_interpolated)
    colorbar
    colormap(jet)
    title(colorbar, 'mGal', 'FontSize', 10);
    title(['column ', num2str(column_num), '-EGM2008'])

    % Save the figure as a PNG file
    saveas(gcf, ['subplot_column_', num2str(column_num), '.png'])

    disp('Mean of Difference GGM and Airborne ')
    mean(Grav_anom - GGM_Gi_interpolated)
end




























