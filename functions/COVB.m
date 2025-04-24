function COVB(KTYPE)
    % ! Determines the type of covariance function to be computed based on KTYPE.
    % ! The covariance function is set as follows:
    % ! - KTYPE = 1: Covariance between the gravity anomaly at P and at Q.
    % ! - KTYPE = 2: Covariance between the longitudinal component of the deflection 
    % !              of the vertical at P and Q.
    % ! - KTYPE = 3: Covariance between the gravity anomaly at P and the height 
    % !              anomaly at Q.
    % ! - KTYPE = 4: Covariance between the longitudinal component of the deflection 
    % !              of the vertical at P and the same at Q.
    % ! - KTYPE = 5: Covariance between the transversal component of the deflection 
    % !              of the vertical at P and the same at Q.
    % ! - KTYPE = 6: Covariance between the longitudinal component of the deflection 
    % !              of the vertical at P and the height anomaly at Q.
    % ! - KTYPE = 7: Covariance between the height anomaly at P and at Q.
    % !
    % ! The value of KTYPE also determines which coefficients (151)-(153) are used
    % ! in the evaluation of the Legendre series and whether differentiation 
    % ! (one or two times) with respect to the variable T is required.
    % ! Logical variables NOTD and NOTDD are used to distinguish between these cases.
    
    % Select covariance type
%     if MODEL
%         return;
%     end
    
    if KTYPE == 1
        IP = 2;
    elseif KTYPE == 2 || KTYPE == 3
        IP = 1;
    else
        IP = 0;
    end
    
    for I = 3:N1
        EPSC(I) = EPS(I) * (((I - 2) * D5 / RBJ) ^ IP);
    end
    NOTD = (KTYPE == 1) || (KTYPE == 3) || (KTYPE == 7);
    NOTDD = (KTYPE ~= 5) && (KTYPE ~= 4);
end