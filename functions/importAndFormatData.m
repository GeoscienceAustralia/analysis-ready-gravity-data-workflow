function [Gravity6D,GravityGradient5D,DEM3D,ZDEM_griddedInterpolant,LongDEMmatrix,LatDEMmatrix,...
    GravityGGM_griddedInterpolant,ZetaGGM_griddedInterpolant,LevellingData3D,...
    ZetaRef_griddedInterpolant,GridRef3D,Coastline,DEM_PARA]=importAndFormatData...
    (GRID_PARA,DEM_PARA,GRAV_PARA,Topo_PARA,COAST_PARA,LEVELLING_PARA,GGM_PARA,GRAV_GRAD_PARA)
    % importAndFormatData import all different data sets needed for LSC.
    %         
    % Input:  GRID_PARA = grid/tiling parameters such as extent (MINLONG,MAXLONG,MINLAT,MAXLAT), buffer, STEP 
    %         DEM_PARA =  DEM data such as filename,num_cols,num_rows
    %         GRAV_PARA =  gravity data such as filename,TypeB ,Grav_Faye_TypeB
    %         COAST_PARA = coast line filename
    %         LEVELLING_PARA = levelling datacomparisons (Lev_eval,filename,Plot_Stat,Compare_To_Existing_Model,Existing_Model, max_diff)
    %         GGM_PARA = GGM file name
    %         GRAV_GRAD_PARA = gravity gradiometry data (filename,TypeB,avail)
    % 
    % Output: Gravity6D = [longitude latitude orthometric_height gravity_anomaly uncertainty(std) flag(1:terrestrial 2:satellite altimetry 3:airborne)]
    %         GravityGradient5D = [longitude latitude height gravity_gradient flag (5:gradiometry)]
    %         DEM3D = [longitude latitude height]
    %         ZDEM_griddedInterpolant = DEM elevation griddedInterpolant
    %         LongDEMmatrix = matrix of longitudes for DEM 
    %         LatDEMmatrix =  matrix of latitudes for DEM 
    %         GravityGGM_griddedInterpolant = GGM gravity griddedInterpolant 
    %         ZetaGGM_griddedInterpolant = GGM height anomaly griddedInterpolant 
    %         LevellingData3D = levelling data [longitude latitude height]
    %         ZetaRef_griddedInterpolant = reference height anomlay griddedInterpolant 
    %         GridRef3D = [LongDEM(:),LatDEM(:),ZDEM(:)]; These are the x,y,z locations the quasigeoid model will be computed at
    %         Coastline= struct with fields lat: long:
    % 
    % Example: see Run_LSC_Drive_RTM
    %
    % Main functions
    % - 
    % Other functions
    % -  importdata
    % -  inpolygon
    %
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11.
    disp('Reading the data sets')
    % Import coastline
    Coastline=importdata(COAST_PARA.filename);
    disp('Gravity')
    Gravity6D=importdata(GRAV_PARA.filename);
    Gravity6D(:,1)=round(Gravity6D(:,1)*60)/60;% make sure we have 1 arc minute data
    Gravity6D(:,2)=round(Gravity6D(:,2)*60)/60;% make sure we have 1 arc minute data
    % this part added for NSW only, for gravity anomaly from gravity gradiometry
    if ~isempty(GRAV_PARA.filename1)
        Gravity5D=importdata(GRAV_PARA.filename1);
        Nanfilter = createNanFilter(Gravity5D(:,4),470,608);% for new NSW file on 22/03/2024
        Gravity5D = Gravity5D(~isnan(Gravity5D(:,4).*Nanfilter), :);
        Gravity5D = [Gravity5D ones(length(Gravity5D), 1) * 3];% add flag column for gravity from gravity gradiometry
        Gravity6D=[Gravity5D;Gravity6D];
    end
    disp('DEM')
    DEM3D=importdata(DEM_PARA.filename);
    % DEM long lat extension: [93 174 -61 -8] 
    LongDEMmatrix=reshape(DEM3D(:,1),DEM_PARA.num_cols,DEM_PARA.num_rows)';
    LatDEMmatrix=reshape(DEM3D(:,2),DEM_PARA.num_cols,DEM_PARA.num_rows)';
    ZDEM=reshape(DEM3D(:,3),DEM_PARA.num_cols,DEM_PARA.num_rows)';
    %disp('Extracting DEM subset') 
    % make sure DEM is bigger than gravity
%     Topo_buffer=Topo_PARA.Rad+GRID_PARA.buffer; 
%     CoordsMM_topo=[GRID_PARA.MINLONG-Topo_buffer,GRID_PARA.MINLAT-Topo_buffer;...
%               GRID_PARA.MINLONG-Topo_buffer,GRID_PARA.MAXLAT+Topo_buffer;...
%               GRID_PARA.MAXLONG+Topo_buffer,GRID_PARA.MAXLAT+Topo_buffer;...
%               GRID_PARA.MAXLONG+Topo_buffer,GRID_PARA.MINLAT-Topo_buffer;...
%               GRID_PARA.MINLONG-Topo_buffer,GRID_PARA.MINLAT-Topo_buffer];
% 
%     DEMin=inpolygon(DEM3D(:,1),DEM3D(:,2),CoordsMM_topo(:,1),CoordsMM_topo(:,2));
%     DEM3D(DEMin==0,:)=[];
    % Computing grid dimensions for one-minute spatial resolution
%     DEM_PARA.num_cols=(max(DEM3D(:,1))-min(DEM3D(:,1)))*60+1;
%     DEM_PARA.num_rows=(max(DEM3D(:,2))-min(DEM3D(:,2)))*60+1;
    % Set the computational grid nodes
%     LongDEMmatrix=reshape(DEM3D(:,1),DEM_PARA.num_cols,DEM_PARA.num_rows)';
%     LatDEMmatrix=reshape(DEM3D(:,2),DEM_PARA.num_cols,DEM_PARA.num_rows)';
%     ZDEM=reshape(DEM3D(:,3),DEM_PARA.num_cols,DEM_PARA.num_rows)';

    ZDEM_griddedInterpolant=griddedInterpolant(LongDEMmatrix',LatDEMmatrix(end:-1:1,:)',ZDEM(end:-1:1,:)');
    
    GridRef3D=[LongDEMmatrix(:),LatDEMmatrix(:),ZDEM(:)];% These are the x,y,z locations the quasigeoid model will be computed at
  
    disp('GGM')
    GGM=importdata(GGM_PARA.filename);
    GravityGGM_griddedInterpolant=griddedInterpolant(GGM.x,GGM.y,GGM.z,GGM.g);
    ZetaGGM_griddedInterpolant=griddedInterpolant(GGM.x,GGM.y,GGM.z,GGM.zeta);
    % Import GPS/leveling
    if LEVELLING_PARA.Compare_To_Existing_Model
        disp('GPS/leveling')
        ZetaRef_griddedInterpolant=importdata(LEVELLING_PARA.Existing_Model);
    else
        ZetaRef_griddedInterpolant=[];
    end
    % Extract data inside computational extents
    disp('Extracting gravity subset') 
    CoordsMM_Grav=[GRID_PARA.MINLONG-GRID_PARA.buffer,GRID_PARA.MINLAT-GRID_PARA.buffer;...
                   GRID_PARA.MINLONG-GRID_PARA.buffer,GRID_PARA.MAXLAT+GRID_PARA.buffer;...
                   GRID_PARA.MAXLONG+GRID_PARA.buffer,GRID_PARA.MAXLAT+GRID_PARA.buffer;...
                   GRID_PARA.MAXLONG+GRID_PARA.buffer,GRID_PARA.MINLAT-GRID_PARA.buffer;...
                   GRID_PARA.MINLONG-GRID_PARA.buffer,GRID_PARA.MINLAT-GRID_PARA.buffer];
            
    Gravin=inpolygon(Gravity6D(:,1),Gravity6D(:,2),CoordsMM_Grav(:,1),CoordsMM_Grav(:,2));
    Gravity6D(Gravin==0,:)=[];
    % import gravity gradiometry data
    if GRAV_GRAD_PARA.avail
    disp('Gravity Gradiometry')
    GravityGradient5D=importdata(GRAV_GRAD_PARA.filename);
%     Grav_grad(:,1)=round(Grav_grad(:,1)*60)/60;% make sure we have 1 arc minute data
%     Grav_grad(:,2)=round(Grav_grad(:,2)*60)/60;% make sure we have 1 arc minute data
%     Grav_grad(:,4)=Grav_grad(:,4);% Change from Eotvos to mgal/m
%     Grav_grad(:,5)=Grav_grad(:,5);% Change from Eotvos to mgal/m
%     this works for Xcalibur_FVD_GDD.mat
      %Nanfilter = createNanFilter(GravityGradient5D(:,4),470,606);
      Nanfilter = createNanFilter(GravityGradient5D(:,4),470,608);%for new NSW file on 22/03/2024
      GravityGradient5D = GravityGradient5D(~isnan(GravityGradient5D(:,4).*Nanfilter), :);
      disp('Extracting gravity gradient subset') 
      GravGradin=inpolygon(GravityGradient5D(:,1),GravityGradient5D(:,2),CoordsMM_Grav(:,1),CoordsMM_Grav(:,2));
      GravityGradient5D(GravGradin==0,:)=[];
    else 
      GravityGradient5D=[];
    end
    % Levelling data
    if LEVELLING_PARA.Lev_eval
        % Buffer these so that we dont include comparisons at the grids edge.
        CoordsMM=[GRID_PARA.MINLONG+GRID_PARA.buffer2,GRID_PARA.MINLAT+GRID_PARA.buffer2;...
                  GRID_PARA.MINLONG+GRID_PARA.buffer2,GRID_PARA.MAXLAT-GRID_PARA.buffer2;...
                  GRID_PARA.MAXLONG-GRID_PARA.buffer2,GRID_PARA.MAXLAT-GRID_PARA.buffer2;...
                  GRID_PARA.MAXLONG-GRID_PARA.buffer2,GRID_PARA.MINLAT+GRID_PARA.buffer2;...
                  GRID_PARA.MINLONG+GRID_PARA.buffer2,GRID_PARA.MINLAT+GRID_PARA.buffer2];
        LevellingData3D=importdata(LEVELLING_PARA.filename);
        Levelling_inpolygon=inpolygon(LevellingData3D(:,1),LevellingData3D(:,2),CoordsMM(:,1),CoordsMM(:,2));
        LevellingData3D=LevellingData3D(Levelling_inpolygon==1,:);
    else
        LevellingData3D=[];
    end
end