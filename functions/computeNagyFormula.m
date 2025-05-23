function prism_out = computeNagyFormula(x,y,z,gravitytype)
    % computeNagyFormula calculates the vertical component of the gravitational potential due to a rectangular 
    % prism at a point outside the prism.
    %
    % Reference: Equation(8) in Nagy, D., Papp, G. & Benedek, J. 
    % The gravitational potential and its derivatives for the prism. 
    % Journal of Geodesy 74, 552–560 (2000). 
    %
    % Input:x = x coordinate vector
    %       y = y coordinate vector
    %       z = z coordinate vector
    %
    % Output: prism_out = vector
    %
    % Example: see computePrismGravity
    %
    % Main functions
    % -  
    % Other functions
    % -   
    % Written by Jack McCubbine
    % Last updated by Neda Darbeheshti
    % Geoscience Australia, 2023-11.

    r = sqrt(x.^2 + y.^2 + z.^2);

    prism_gradient = atan((x.*y)./(z.*r)); 

    switch gravitytype
        case 'gg' 
            prism_out = prism_gradient;  
        case 'g'
            prism_out = -x.*log(y+r) - y.*log(x+r) + z.*prism_gradient;  
        otherwise
            warning('Unexpected gravity type. No prism created.')
    end

end