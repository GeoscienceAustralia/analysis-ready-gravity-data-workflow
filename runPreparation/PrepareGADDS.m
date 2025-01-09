clear
% Height is ground elevation in GADDS
australia2016GADDS = readtable('Data/GRAVITY/GADDS/gravity_may2016.dat');
VicNSWNov2024GADDS = readtable('Data/GRAVITY/GADDS/NSW_Vic_Nov2024.csv');


% Select the 9 columns you want 
selectedColumns = australia2016GADDS(:, [4, 5, 6, 7, 8, 9, 18, 19,str2double(:,22)]);

for i = 1:width(selectedColumns)
    selectedColumns{:, i} = str2double(selectedColumns{:, i});
end




% Divide the last and middle columns by 10
selectedColumns{:, 9} = selectedColumns{:, 9} / 10;
selectedColumns{:, 6} = selectedColumns{:, 6} / 10;







