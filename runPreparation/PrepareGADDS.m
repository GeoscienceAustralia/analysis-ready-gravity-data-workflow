clear
%% Import the raw GADDS gravity data
GADDS=fopen('Data/GRAVITY/GADDS/gravity_may2016.dat');
TEMP=strsplit(fgetl(GADDS));

t=0;
%tic
data=zeros(1800000,9);
for k=1:length(data(:,1))
    try
                   %[Longitude, Latitude, Height, height, N, Gravity, Horizontal Uncertainty, Vertical Uncertainty]
        data(k,:)=[str2double(TEMP{5}),str2double(TEMP{6}),str2double(TEMP{7}),str2double(TEMP{8}),str2double(TEMP{9}),str2double(TEMP{10})/10,str2double(TEMP{19}),str2double(TEMP{20}),str2double(TEMP{23})/10];
        TEMP=strsplit(fgetl(GADDS));
        t=k;
        if mod(k,1000)==0
            toc
            tic
            disp(k)
        end
    catch
    end
end
data(t+1:end,:)=[];
%toc
fclose(GADDS);


australia2016GADDS = readtable('Data/GRAVITY/GADDS/gravity_may2016.dat');
VicNSWNov2024GADDS = readtable('Data/GRAVITY/GADDS/NSW_Vic_Nov2024.csv');
