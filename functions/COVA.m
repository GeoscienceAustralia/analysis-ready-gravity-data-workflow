function COVA(EPS, N1)
    %     ! This subroutine computes one of seven different covariances (see below)
    %     ! using the anomaly degree-variance model given through the values of 
    %     ! Table Seven and Equation (68). The quantity S in the table is 
    %     ! here called SE.
    %     
    %     ! There are three entries to the subroutine, which have to be called in
    %     ! the sequence COVA, COVB, and COVC.
    %     
    %     ! By calling COVA, the kind of covariance function to be used is determined.
    %     ! There are three possibilities:
    %     ! (1) The covariance model four (Equations (130)-(132) and (136)-(139))
    %     !     is used without modifications. In this case, EPS will be a dummy
    %     !     array and N1 must be equal to one. The logical variable MODEL will
    %     !     get the value TRUE in this case.
    %     ! (2) A number (N1) of the anomaly degree-variances (degree zero to N1-1)
    %     !     are put equal to empirically determined degree-variances. The degree-
    %     !     variance of degree K will have to be stored in EPS(K+1) (in units
    %     !     of MGAL**2).
    %     ! (3) The degree-variances of degree zero to N=N1-1 are put equal to zero,
    %     !     and the others are the same as described above. This means that an
    %     !     N'th order local covariance function will be computed. In this case,
    %     !     EPS must have NZ zero values stored.
    %     
    %     ! In all cases, N1 must be less than 300 and EPS must have dimension N1.
    % Constants
    RE = 6371.0e3;
    GM = 3.98e14;
    A = 425.28;    % from Table Seven Tscherning(1974)
    SE = 0.999617; % from Table Seven Tscherning(1974)
    B = 24.0;      % from Table Seven Tscherning(1974)
    IB1 = 25;
    IB2 = 26;
    IBM1 = 23;
    EPSC = zeros(1, 300); % Assuming EPSC is a vector of size 300
    EPSC(1) =[] ;
    EPSC(2) = [];
    D0 = 3 * 0.0;
    D1 = 1.0;
    D2 = 2.0;
    D3 = 3.0;
    D4 = 4.0;
    D5 = 1.0e5;
    RADSEC = 206264.806; 

    IB12 = IB1 * IB2;
    RADSE2 = RADSEC^2;
    RE2 = RE^2;
    RBJ2 = RE2 * SE;
    RBJ = sqrt(RBJ2);
    AM = A / D5;
    AM2 = AM / D5;
 %  A IS IN UNITS OF MGAL**2, AM IN UNITS OF MGAL*M/SEC AND AM2 IN 
 %  UNITS OF (M/SEC)**2. RBJ IS THE RADIUS OF THE BJERHAMMAK-SPHERE. 
    MODEL = (N1 == 1);
    if MODEL
        return;
    end
  % WE WILL NOW COMPUTE THE MODIFIED (POTENTIAL) DEGREE-VARIANCES, CF. 
  % EQUATION (151).
    if N1 < 3
        return;
    end
    
    for I = 3:N1
        RI = double(I - 1);
        if I == 3
            EPS(3) = EPS(3) * RBJ2 * 1.0e-10;
        end
        if I > 3
            EPS(I) = RBJ2 * (EPS(I) / ((RI - D1)^2)) * 1.0e-10 - AM2 / ((RI - D1) * (RI - D2) * (RI + 8));
        end
    end
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
    if MODEL
        return;
    end
    
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