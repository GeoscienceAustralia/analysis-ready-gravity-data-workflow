function comparetoLevellingData(LEVELLING_PARA,GRID_PARA,LongDEM,LatDEM,Geoid,Grid_res_geoid,Lev, ...
    REFERENCE_GEOID_Zetai,GGM_Zeta_griddedInterpolant,Weights,Coastline)
    % This function compares gravimetric geoid to geometric geoid.
    %
    % Input:  sphericalDistance = vector of spherical distances
    %         empiricalCovariance = vector of empirical covariance
    %         maxOrder = maximum order for Legendre polynomials
    %         minOrder = minimum order for Legendre polynomials
    % 
    % Output: 
    %           
    % Example: see computeLSC
    %
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11. 

    [Longm_i,Latm_i]=meshgrid(unique(LongDEM(1,:)),flip(unique(LatDEM(:,2))));

    Geoid_LSC_Function=griddedInterpolant(Longm_i',Latm_i(end:-1:1,:)',Geoid(end:-1:1,:)');

    Vals_Lev=Geoid_LSC_Function(Lev(:,1),Lev(:,2))-Lev(:,3);    
    Vals_GGM=GGM_Zeta_griddedInterpolant(Lev(:,1),-Lev(:,2),Lev(:,2)*0)-Lev(:,3);

    Vals_Levog=Vals_Lev;
    Vals_GGMog=Vals_GGM;

    Vals_Levog(abs(Vals_Lev-median(Vals_Lev,"omitnan"))>LEVELLING_PARA.max_diff)=[];
    Vals_GGMog(abs(Vals_Lev-median(Vals_Lev,"omitnan"))>LEVELLING_PARA.max_diff)=[];

    disp('Levelling GGM')
    disp('mean min max std')
    disp([mean(Vals_GGM(isnan(Vals_Levog)==0),"omitnan"),min(Vals_GGM(isnan(Vals_Levog)==0)),max(Vals_GGM(isnan(Vals_Levog)==0)),std(Vals_GGM(isnan(Vals_Levog)==0),"omitnan")])
    disp('Levelling LSC')
    disp('mean min max std')
    disp([mean(Vals_Levog,"omitnan"),min(Vals_Levog),max(Vals_Levog),std(Vals_Levog,"omitnan")])

    if LEVELLING_PARA.Compare_To_Existing_Model
        Vals_REFERENCE_GEOID=REFERENCE_GEOID_Zetai(Lev(:,1),Lev(:,2))-Lev(:,3);
        Vals_REFERENCE_GEOIDog=Vals_REFERENCE_GEOID;
        Vals_REFERENCE_GEOIDog(abs(Vals_Lev-median(Vals_Lev,"omitnan"))>0.15)=[];
        Vals_REFERENCE_GEOIDog(isnan(Vals_Levog))=NaN;
        disp('Levelling Reference Geoid')
        disp('mean min max std')
        disp([mean(Vals_REFERENCE_GEOIDog,"omitnan"),min(Vals_REFERENCE_GEOIDog),max(Vals_REFERENCE_GEOIDog),std(Vals_REFERENCE_GEOIDog,"omitnan")])
    end
   % plot levelling data
    if LEVELLING_PARA.Plot_Stats
        plotLevellingData(LongDEM,LatDEM,Lev,Vals_Levog,Vals_REFERENCE_GEOID,Vals_REFERENCE_GEOIDog,Vals_GGM,Vals_Lev, ...
        Grid_res_geoid,Weights,GGM_Zeta_griddedInterpolant,REFERENCE_GEOID_Zetai,Coastline,GRID_PARA,LEVELLING_PARA)    
    end
end