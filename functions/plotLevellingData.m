function plotLevellingData(LongDEM,LatDEM,Lev,Vals_Levog,Vals_REFERENCE_GEOID,Vals_REFERENCE_GEOIDog,Vals_GGM,Vals_Lev, ...
    Grid_res_geoid,Weights,GGM_Zeta_griddedInterpolant,REFERENCE_GEOID_Zetai,Coastline,GRID_PARA,LEVELLING_PARA) 

            figure(7)
            clf
            histogram(Vals_Levog(isnan(Vals_Levog)==0)-median(Vals_Levog(isnan(Vals_Levog)==0),"omitnan"),20)

            figure(8)
            clf
            hold on
            plot(sort(Vals_Levog(isnan(Vals_Levog)==0))-median(Vals_Levog(isnan(Vals_Levog)==0),"omitnan"),0:1/length(Vals_Levog(isnan(Vals_Levog)==0)):1-1/length(Vals_Levog(isnan(Vals_Levog)==0)),'g')
            if LEVELLING_PARA.Compare_To_Existing_Model
            plot(sort(Vals_REFERENCE_GEOIDog(isnan(Vals_Levog)==0))-median(Vals_REFERENCE_GEOIDog(isnan(Vals_Levog)==0),"omitnan"),0:1/length(Vals_Levog(isnan(Vals_Levog)==0)):1-1/length(Vals_Levog(isnan(Vals_Levog)==0)),'r--')
            legend('LSC','Reference Geoid model')
            else
            legend('LSC')
            end
            title('Cumulative distribution plot')
            grid on

            figure(9)
            clf
            hold on
            plot(sort(Vals_GGM(isnan(Vals_Levog)==0))-median(Vals_GGM(isnan(Vals_Levog)==0),"omitnan"),0:1/length(Vals_GGM(isnan(Vals_Levog)==0)):1-1/length(Vals_GGM(isnan(Vals_Levog)==0)),'k')
            plot(sort(Vals_Levog(isnan(Vals_Levog)==0))-median(Vals_Levog(isnan(Vals_Levog)==0),"omitnan"),0:1/length(Vals_Levog(isnan(Vals_Levog)==0)):1-1/length(Vals_Levog(isnan(Vals_Levog)==0)),'g')
            if LEVELLING_PARA.Compare_To_Existing_Model
            plot(sort(Vals_REFERENCE_GEOIDog(isnan(Vals_Levog)==0))-median(Vals_REFERENCE_GEOIDog(isnan(Vals_Levog)==0),"omitnan"),0:1/length(Vals_Levog(isnan(Vals_Levog)==0)):1-1/length(Vals_Levog(isnan(Vals_Levog)==0)),'r--')
            legend('GGM','LSC','Reference Geoid model')
            else
            legend('GGM','LSC')
            end
            title('Cumulative distribution plot')
            grid on

            figure(10)
            clf
            hold on
            for jj=1:9
            plot(Coastline.long{jj},Coastline.lat{jj},'k')
            end
            scatter(Lev(:,1),Lev(:,2),50,Vals_Lev-median(Vals_Lev,"omitnan"),'filled')
            colorbar
            colormap(jet)
            axis([GRID_PARA.MINLONG GRID_PARA.MAXLONG GRID_PARA.MINLAT GRID_PARA.MAXLAT])
            ax=gca;
            ax.PlotBoxAspectRatio =[1 abs(cos(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT])*pi/180)) 1];
            caxis([-0.25 0.25])
            title('LSC Geoid minus GPS-levelling data')
            drawnow

            figure(11)
            clf
            hold on
            for jj=1:9
                plot(Coastline.long{jj},Coastline.lat{jj},'k')
            end
            scatter(Lev((isnan(Vals_Lev)==0),1),Lev((isnan(Vals_Lev)==0),2),50,Vals_GGM(isnan(Vals_Lev)==0)-median(Vals_GGM((isnan(Vals_Lev)==0)),"omitnan"),'filled')
            colorbar
            colormap(jet)
            axis([GRID_PARA.MINLONG GRID_PARA.MAXLONG GRID_PARA.MINLAT GRID_PARA.MAXLAT])
            ax=gca;
            ax.PlotBoxAspectRatio =[1 abs(cos(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT])*pi/180)) 1];
            caxis([-0.25 0.25])
            title('GGM minus GPS-levelling data')
            drawnow

            if LEVELLING_PARA.Compare_To_Existing_Model
                figure(12)
                clf
                hold on
                for jj=1:9
                    plot(Coastline.long{jj},Coastline.lat{jj},'k')
                end
                scatter(Lev((isnan(Vals_Lev)==0),1),Lev((isnan(Vals_Lev)==0),2),50,Vals_REFERENCE_GEOID(isnan(Vals_Lev)==0)-median(Vals_REFERENCE_GEOID((isnan(Vals_Lev)==0)),"omitnan"),'filled')
                colorbar
                colormap(jet)
                axis([GRID_PARA.MINLONG GRID_PARA.MAXLONG GRID_PARA.MINLAT GRID_PARA.MAXLAT])
                ax=gca;
                ax.PlotBoxAspectRatio =[1 abs(cos(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT])*pi/180)) 1];
                caxis([-0.25 0.25])
                title('Reference Geoid minus GPS-levelling data')
                drawnow

                figure(13)
                clf
                hold on
                imagesc(LongDEM(1,:),LatDEM(:,1),(REFERENCE_GEOID_Zetai(LongDEM,LatDEM)-GGM_Zeta_griddedInterpolant(LongDEM,-LatDEM,LatDEM*0))-mean(mean(REFERENCE_GEOID_Zetai(LongDEM,LatDEM)-GGM_Zeta_griddedInterpolant(LongDEM,-LatDEM,LatDEM*0))))
                for jj=1:9
                    plot(Coastline.long{jj},Coastline.lat{jj},'k')
                end
                colorbar
                colormap(jet)
                axis([GRID_PARA.MINLONG GRID_PARA.MAXLONG GRID_PARA.MINLAT GRID_PARA.MAXLAT])
                ax=gca;
                ax.PlotBoxAspectRatio =[1 abs(cos(mean([GRID_PARA.MINLAT GRID_PARA.MAXLAT])*pi/180)) 1];
                title('Residual Reference Geoid')
                caxis([-2*std(Grid_res_geoid(:)./Weights(:),"omitnan") 2*std(Grid_res_geoid(:)./Weights(:),"omitnan")])
                drawnow

            end