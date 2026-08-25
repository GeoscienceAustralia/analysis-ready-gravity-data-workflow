function plotGGMfunctionals(GRID_PARA,GRAV_GRAD_PARA,OUTPUT_PARA,Grav,Grav_grad,LongDEM,LatDEM, ...
    GGM_Gravity_griddedInterpolant,GGM_Zeta_griddedInterpolant,ZDEM_griddedInterpolant,Coastline)

    GGM_Gravity = GGM_Gravity_griddedInterpolant(Grav(:,1),-Grav(:,2),Grav(:,3)-ZDEM_griddedInterpolant(Grav(:,1),Grav(:,2)));
    
    Grav(:,4) = Grav(:,4) - GGM_Gravity;
    
    plotCustomScatter(Grav(:,1),Grav(:,2),GGM_Gravity,GRID_PARA,'GGM gravity','mGal',Coastline,[],OUTPUT_PARA.plotsFolder)
        
   if GRAV_GRAD_PARA.avail

       GGM_GravityGradient=(GGM_Gravity_griddedInterpolant(Grav_grad(:,1),-Grav_grad(:,2), ...
                   Grav_grad(:,3)-ZDEM_griddedInterpolant(Grav_grad(:,1),Grav_grad(:,2))-0.5)-...
                   GGM_Gravity_griddedInterpolant(Grav_grad(:,1),-Grav_grad(:,2), ...
                   Grav_grad(:,3)-ZDEM_griddedInterpolant(Grav_grad(:,1),Grav_grad(:,2))+0.5));
        
       Grav_grad(:,4) = Grav_grad(:,4) - GGM_GravityGradient; 
    
       plotCustomScatter(Grav_grad(:,1),Grav_grad(:,2),GGM_GravityGradient,GRID_PARA,'GGM gravity gradient','mGal/m',Coastline,[],OUTPUT_PARA.plotsFolder)
          
    end

% common variables for plotting
axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;
figure('Name','GGM','NumberTitle','off'); 
clf
hold on
imagesc(LongDEM(1,:),LatDEM(:,1),GGM_Zeta_griddedInterpolant(LongDEM,-LatDEM,LatDEM*0))
customizeMap('GGM Zeta','m',Coastline,axisLimits)
%caxis([-0.1 0.1])
saveas(gcf,[OUTPUT_PARA.plotsFolder,'GGMzeta','.png']) 

end


% EGM_PARA.filename='Data/GGM/EGM2008_For_Gridded_Int.mat';
% 
% EGM=importdata(EGM_PARA.filename);
% 
% EGM_Zeta_griddedInterpolant=griddedInterpolant(EGM.x,EGM.y,EGM.z,EGM.zeta);
% 
% EGM_ZetaonDEM=EGM_Zeta_griddedInterpolant(LongDEM,-LatDEM,LatDEM*0);
% 
% ZDegEGM=mean(mean(REFERENCE_Zeta_griddedInterpolant(LongDEM,LatDEM)-EGM_ZetaonDEM));
% 
% EGMresAGQG=REFERENCE_Zeta_griddedInterpolant(LongDEM,LatDEM)-EGM_ZetaonDEM;
% 
% ZDegEGMLSCAGQG=mean(mean(Grid_res_geoid_w-EGM_ZetaonDEM));
% 
% EGMresLSCAGQG=Grid_res_geoid_w-EGM_ZetaonDEM;
% 
% % common variables for plotting
% axisLimits.latMeanCosine=abs(cos(deg2rad(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT]))));
% axisLimits.lonMinLimit=GRID_PARA.MINLONG-GRID_PARA.buffer;
% axisLimits.lonMaxLimit=GRID_PARA.MAXLONG+GRID_PARA.buffer;
% axisLimits.latMinLimit=GRID_PARA.MINLAT-GRID_PARA.buffer;
% axisLimits.latMaxLimit=GRID_PARA.MAXLAT+GRID_PARA.buffer;
% 
% % plot residualGeoidvsAGQG
% figure('Name','MosaicTiles','NumberTitle','off'); 
% clf
% subplot(1,2,1);
% hold on
% imagesc(LongDEM(1,:),LatDEM(:,1),EGMresLSCAGQG-ZDegEGMLSCAGQG)
% customizeMap('ResEGM2008 LSC AGQG','m',Coastline,axisLimits)
% caxis([-0.1 0.1])
%  
% subplot(1,2,2);
% %figure
% hold on
% imagesc(LongDEM(1,:),LatDEM(:,1),EGMresAGQG-ZDegEGM)
% customizeMap('ResEGM2008 2022 AGQG','m',Coastline,axisLimits) 
% caxis([-0.1 0.1])
% saveas(gcf,[OUTPUT_PARA.plotsFolder,'MosaicTiles','EGM2008residualLSCvs2022AGQG','.png']) 
% 
% diary off