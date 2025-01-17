close all
clear
% Height is ground elevation in GADDS

% Read the CSV file into a matrix (numeric values as double)
NSW_Vic_Data = readmatrix('Data/GRAVITY/GADDS/NSW_Vic_Nov2024.csv');

%[Longitude, Latitude, Height, height, N, Gravity, Horizontal Uncertainty, Vertical Uncertainty]
selectedColumns = NSW_Vic_Data(:, [6, 7, 14, 22, 49, 29,11,19]);























