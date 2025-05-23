function normalGravity = computeNormalGravity (lat)

    % compute SOMIGLIANA CLOSED NORMAL GRAVITY (MGAL)
    % reference: Eqa(4) in W. E. FEATHERSTONE and M.C. DENTITH (1997)
    % Input:  lat = latitude in radians
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2024-11.

    constants                                       % load constants

    normalGravity = AbsoluteGravityEquator_mgal*(1+NormalGravityConstant*(sin(lat).^2) ...
            )./sqrt(1-EarthEccentricitySquared*(sin(lat).^2));
end

