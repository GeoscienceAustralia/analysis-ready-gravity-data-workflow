function atmospheric_corr = computeAtmosphericGravityCorrection (ht)
    %  COMPUTES THE ATMOSPHERIC GRAVITY CORRECTION (MGAL)
    % Input:  ht =  height
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2024-11.

    atmospheric_corr = 0.871 - 1.0298*(10^-4)*ht + 5.3105*(10^-9)*(ht.^2) - 2.1642*(10^-13)*(ht.^3);
end



