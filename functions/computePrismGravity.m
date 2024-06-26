function [terrc] = computePrismGravity(x1,x2,y1,y2,z1,z2,rho,gravitytype)
    %
    % computePrismGravity performs terrain contribution to the gravity using the Nagy 
    % prism formula. 
    % To calculate the terrain contribution to the gravity for a prism mass
    % with coordinates in the x,y,z plane and
    %.........................................
    %............______________(x2,y2,z2).....
    %.........../!............/!..............
    %........../_!___________/.!..............
    %..........!.!..........!..!..............
    %..........!.!..........!..!..............  
    %..........!.!..........!..!..............  
    %..........!.!..........!..!..............  
    %..........!.!(x1,y2,z1)!..!(x2,y2,z1)....  
    %..........!./..........!./............... 
    %(x1,y1,z1)!/___________!/(x2,y1,z1)......
    % 
    % with density rho;
    % the terrain gravity dT/dz is given by
    % dT/dz = computePrismGravity(x1,x2,y1,y2,z1,z2,rho)
    % 
    % The function can handle vectors of coordinates so that correction
    % for mutiple pieces of terrain can also be calculated
    % LHS=(xp-x1) RHS=(xp-x2)
    % Front=(yp-y1) back=(yp-y2)
    % Bottom=(zp-z1) Top=(zp-z2)
    %
    % Input:    x1= vector
    %           x2= vector
    %           y1= vector
    %           y2= vector
    %           z1= number
    %           z2= vector
    %           rho= number, assumed density in g/cm^3.
    % 
    % Output:   terrc= vector 
    %           
    % Example: see computeTerrainCorrection
    %
    % Main functions
    % - computeNagyFormula
    % - 
    % Other functions
    % -   
    %
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11.
  
    constants                                       % load constants
    % gravitational constant in dyn⋅cm2⋅g−2
    %bigG = 6.6720e-08; 
    % convert bigG to mGal
    bigG_mgal=bigG*(10^5);
    %
    terrc=((computeNagyFormula(x2,y2,z2,gravitytype)-computeNagyFormula(x2,y2,z1,gravitytype))-(computeNagyFormula(x2,y1,z2,gravitytype)-computeNagyFormula(x2,y1,z1,gravitytype)))-...
          ((computeNagyFormula(x1,y2,z2,gravitytype)-computeNagyFormula(x1,y2,z1,gravitytype))-(computeNagyFormula(x1,y1,z2,gravitytype)-computeNagyFormula(x1,y1,z1,gravitytype)));
    %   
    terrc=bigG_mgal.*rho.*terrc;
end

