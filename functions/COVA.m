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

    
end