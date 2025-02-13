
% Data matrix
data = [
    -52.833, 8, 0, 0.000, 93, 2, 0.000, -1.34, -11.73;
    -52.726, 8, 0, 0.000, 93, 3, 0.000, -1.12, -12.12;
    -52.617, 8, 0, 0.000, 93, 4, 0.000, -1.01, -12.24;
    -52.508, 8, 0, 0.000, 93, 5, 0.000, -0.90, -12.24;
    -52.399, 8, 0, 0.000, 93, 6, 0.000, -0.90, -12.18;
    -52.291, 8, 0, 0.000, 93, 7, 0.000, -1.01, -12.07;
    -52.184, 8, 0, 0.000, 93, 8, 0.000, -1.12, -11.79;
    -52.081, 8, 0, 0.000, 93, 9, 0.000, -1.12, -11.28;
    -51.983, 8, 0, 0.000, 93, 10, 0.000, -1.12, -10.72;
    -51.890, 8, 0, 0.000, 93, 11, 0.000, -1.23, -9.99;
    -51.805, 8, 0, 0.000, 93, 12, 0.000, -1.12, -9.21;
    -51.726, 8, 0, 0.000, 93, 13, 0.000, -1.23, -8.53;
    -51.653, 8, 0, 0.000, 93, 14, 0.000, -1.45, -7.80;
    -51.587, 8, 0, 0.000, 93, 15, 0.000, -1.45, -7.07;
    -51.527, 8, 0, 0.000, 93, 16, 0.000, -1.57, -6.34;
    -51.474, 8, 0, 0.000, 93, 17, 0.000, -1.57, -5.61;
    -51.427, 8, 0, 0.000, 93, 18, 0.000, -1.57, -5.00;
    -51.385, 8, 0, 0.000, 93, 19, 0.000, -1.45, -4.55;
    -51.346, 8, 0, 0.000, 93, 20, 0.000, -1.45, -4.27;
    -51.309, 8, 0, 0.000, 93, 21, 0.000, -1.45, -4.21;
    -51.271, 8, 0, 0.000, 93, 22, 0.000, -1.57, -4.27;
    -51.233, 8, 0, 0.000, 93, 23, 0.000, -1.57, -4.38;
    -51.193, 8, 0, 0.000, 93, 24, 0.000, -1.68, -4.66;
    -51.150, 8, 0, 0.000, 93, 25, 0.000, -1.68, -5.11;
    -51.102, 8, 0, 0.000, 93, 26, 0.000, -1.79, -5.56;
    -51.051, 8, 0, 0.000, 93, 27, 0.000, -1.79, -5.89;
    -50.997, 8, 0, 0.000, 93, 28, 0.000, -1.90, -6.06;
    -50.943, 8, 0, 0.000, 93, 29, 0.000, -2.01, -5.95;
    -50.891, 8, 0, 0.000, 93, 30, 0.000, -2.01, -5.67;
    -50.842, 8, 0, 0.000, 93, 31, 0.000, -2.01, -5.39;
    -50.795, 8, 0, 0.000, 93, 32, 0.000, -2.01, -5.22;
    -50.749, 8, 0, 0.000, 93, 33, 0.000, -2.01, -5.05;
    -50.705, 8, 0, 0.000, 93, 34, 0.000, -2.01, -4.94;
    -50.661, 8, 0, 0.000, 93, 35, 0.000, -1.90, -4.94;
    -50.617, 8, 0, 0.000, 93, 36, 0.000, -1.90, -4.94
];


% Define the data
latitudes = [-53.031, -52.935, -52.833, -52.726, -52.617, -52.508, -52.399, -52.291, -52.184, -52.081, -51.983, -51.890, -51.805, -51.726, -51.653, -51.587, -51.527, -51.474, -51.427, -51.385, -51.346, -51.309, -51.271, -51.233, -51.193, -51.150];
longitudes = 93 * ones(size(latitudes));
values1 = [-1.68, -1.45, -1.34, -1.12, -1.01, -0.90, -0.90, -1.01, -1.12, -1.12, -1.12, -1.23, -1.12, -1.23, -1.45, -1.45, -1.57, -1.57, -1.57, -1.45, -1.45, -1.45, -1.57, -1.57, -1.68, -1.68];
values2 = [-10.78, -11.11, -11.73, -12.12, -12.24, -12.24, -12.18, -12.07, -11.79, -11.28, -10.72, -9.99, -9.21, -8.53, -7.80, -7.07, -6.34, -5.61, -5.00, -4.55, -4.27, -4.21, -4.27, -4.38, -4.66, -5.11];

addpath('outputs');

% Open the file for writing
fileID = fopen('outputs/output.txt', 'w');


% Loop through data and write to the file without commas
for i = 1:size(data, 1)
    fprintf(fileID, 'GEO  %.3f  S  %2d  %2d  %.3f  E  %3d  %2d  %.3f  %.2f  %.2f\n', ...
            data(i, 1), data(i, 2), data(i, 3), data(i, 4), data(i, 5), data(i, 6), data(i, 7), data(i, 8), data(i, 9));
%           value1,     lat_deg,    lat_min,    lat_sec,    lon_deg,    lon_min,    lon_sec,    value2,     value3
end
 
% Close the file
fclose(fileID);







