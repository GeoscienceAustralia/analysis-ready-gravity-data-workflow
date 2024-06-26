function [gravity_correction]=computeTerrainCorrection(ZDEM_topo,LatDEM_topo,LongDEM_topo,Lat_CP,Long_CP,H_CP,Rad,rho,gravityType)
    % This function computes the terrain correction. 
    %
    % Input:    ZDEM_topo = elevation matrix
    %           LatDEM_topo = latitude matrix
    %           LongDEM_topo = longitude matrix
    %           Lat_CP = latitude vector for computation point
    %           Long_CP = longitude vector for computation point
    %           H_CP =  height vector for computation point
    %           Rad = number, Radius out to which to compute the effects in degrees.
    %           rho = number, assumed density in g/cm^3.
    %           gravityType = character, 'g' for gravity
    %                                   'gg' for gravity gradient
    %         
    % Output:   gravity_correction= vector of gravity correction
    %           
    % Example: see computeTerrainEffect
    %
    % Main functions
    % - computePrismGravity
    % - plotTerrainCorrection
    % Other functions
    % -   
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11.
    constants                          % load constants
    % define the differential element for Nagy Prism formula
    min2deg=1/60;
    dy=(deg2meter*min2deg)/2;          % half a minute differential element
    % Topography gridded interpolant function
    ZDEM_topo(ZDEM_topo<0)=0;
    H_CP(H_CP<0)=0;
    % Compute residual topographic effect and potential    
    % Initialise
    gravity_correction=H_CP*NaN;
    IND=1:length(Long_CP(:));
    %             
    for longi=min(LongDEM_topo(:))+Rad:2*Rad:max(LongDEM_topo(:))
        
        for lati=min(LatDEM_topo(:))+Rad:2*Rad:max(LatDEM_topo(:))
      
            LatDEM_topo_loc=LatDEM_topo(LatDEM_topo>lati-2*Rad & LatDEM_topo<lati+2*Rad & LongDEM_topo>longi-2*Rad & LongDEM_topo<longi+2*Rad);
            LongDEM_topo_loc=LongDEM_topo(LatDEM_topo>lati-2*Rad & LatDEM_topo<lati+2*Rad & LongDEM_topo>longi-2*Rad & LongDEM_topo<longi+2*Rad);
            ZDEM_topo_loc=ZDEM_topo(LatDEM_topo>lati-2*Rad & LatDEM_topo<lati+2*Rad & LongDEM_topo>longi-2*Rad & LongDEM_topo<longi+2*Rad);
            
            Lat_CP_loc=Lat_CP(Lat_CP>=lati-Rad & Lat_CP<=lati+Rad & Long_CP>=longi-Rad & Long_CP<=longi+Rad);
            Long_CP_loc=Long_CP(Lat_CP>=lati-Rad & Lat_CP<=lati+Rad & Long_CP>=longi-Rad & Long_CP<=longi+Rad);
            H_CP_loc=H_CP(Lat_CP>=lati-Rad & Lat_CP<=lati+Rad & Long_CP>=longi-Rad & Long_CP<=longi+Rad);
            IND_loc=IND(Lat_CP>=lati-Rad & Lat_CP<=lati+Rad & Long_CP>=longi-Rad & Long_CP<=longi+Rad);
            
            % Convert degrees to radians
            LongDEM_topo_locRadian = deg2rad (LongDEM_topo_loc);
            Long_CP_locRadian = deg2rad (Long_CP_loc);
            LatDEM_topo_locRadian=deg2rad (LatDEM_topo_loc);
            Lat_CP_locRadian=deg2rad (Lat_CP_loc);
             
            for k=1:length(H_CP_loc(:))
    
                 %distanceFromCP=sqrt((LongDEM_topo_loc-Long_CP_loc(k)).^2+(LatDEM_topo_loc-Lat_CP_loc(k)).^2);
                 % should be tested
                 distanceFromCP = haversine(LatDEM_topo_locRadian(:), LongDEM_topo_locRadian(:), Lat_CP_locRadian(k), Long_CP_locRadian(k));
                 distanceFromCP=rad2deg(distanceFromCP);
                 
                 ID=1:length(LongDEM_topo_loc(:));
                 IDr=ID(distanceFromCP(:)<=Rad);
                 H_rmt=ZDEM_topo_loc(IDr);
            
                  if sum(H_rmt)~=0
                     Long_rmt=LongDEM_topo_loc(IDr);
                     Lat_rmt=LatDEM_topo_loc(IDr);
    
                     h_top=H_rmt;   
                     H_CP_P=H_CP_loc(k);
    
                     x=(Long_rmt-Long_CP_loc(k))*deg2meter*cos(deg2rad(Lat_CP_loc(k)));
                     y=(Lat_rmt-Lat_CP_loc(k))*deg2meter;
    
                     dx=(deg2meter*cos(deg2rad(Lat_CP_loc(k)))*min2deg)/2;
    
                     gravity_correction(IND_loc(k))=sum(-computePrismGravity(x-dx,x+dx,y-dy,y+dy,-H_CP_P,h_top-H_CP_P,rho,gravityType));
                
                  end
            
            end
    
        end
        
    end
            gravity_correction(isnan(gravity_correction))=0;
 end
