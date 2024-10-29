function fac = computeSecondOrderFreeAirCorrection(lat,ht,gamma)
    % 2nd order Free-Air Correction 
    % reference: GRAV-D General Airborne Gravity Data User Manual, page 23
    % Input:  lat = latitude in radians
    %         ht =  height
    %         gamma = Normal gravity, the absolute gravity of the Earth
    %         approximating ellipsoid(Moritz 1980)
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2024-11.

    constants                                       % load constants
    
    EarthMajorAxis = EarthMajorAxis*10^3;            % km to m
    
    sphi2 = sin(lat).^2;
    
    c1 = ((2.*gamma)./EarthMajorAxis);
    
    term1 = (1+flattening+GravToCentrifugalRatio_Equator-2*flattening*sphi2).*ht;
    
    c2 = (3.*AbsoluteGravityEquator_mgal./(EarthMajorAxis.^2));

    term2 =ht.^2;
    
    fac = (c1.*term1)-(c2.*term2);
end









