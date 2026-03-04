clear
close all
outputName='NQueen1deg';%'NQueenfrom16000'
OUTPUT_PARA.Grids_name=['outputs/Grids',outputName,'/'];
OUTPUT_PARA.Tiles_dir_name=['outputs/ResidualTiles',outputName,'/'];
OUTPUT_PARA.plotsFolder=['outputs/Grids',outputName,'/',date,outputName];


Files=dir(OUTPUT_PARA.Tiles_dir_name);
Files(1:2)=[];


for k=1:length(Files)

    Tile_Data=importdata([OUTPUT_PARA.Tiles_dir_name,Files(k).name]);

    COV_PARA_RTMA(k)=Tile_Data.COV_PARA_RTM.A;
    COV_PARA_RTMB(k)=Tile_Data.COV_PARA_RTM.B;
    COV_PARA_FayeA(k)=Tile_Data.COV_PARA_Faye.A;
    COV_PARA_FayeB(k)=Tile_Data.COV_PARA_Faye.B;
    
end

figure
plot (COV_PARA_FayeA,'.')
title('COVPARAFayeA')
saveas(gcf,[OUTPUT_PARA.plotsFolder,'COVPARAFayeA','.png']) 
figure
plot (COV_PARA_FayeB,'.')
title('COVPARAFayeB')
saveas(gcf,[OUTPUT_PARA.plotsFolder,'COVPARAFayeB','.png']) 
figure
plot (COV_PARA_RTMA,'.')
title('COVPARARTMA')
saveas(gcf,[OUTPUT_PARA.plotsFolder,'COVPARARTMA','.png']) 
figure
plot (COV_PARA_RTMB,'.')
title('COVPARARTMB')
saveas(gcf,[OUTPUT_PARA.plotsFolder,'COVPARARTMB','.png']) 

