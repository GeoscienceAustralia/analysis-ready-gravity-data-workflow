function COV = COVA(EPS, N1)
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
    %     ! (2) A number (NI) of the anomaly degree-variances (degree zero to NI-1)
    %     !     are put equal to empirically determined degree-variances. The degree-
    %     !     variance of degree K will have to be stored in EPS(K+1) (in units
    %     !     of MGAL**2).
    %     ! (3) The degree-variances of degree zero to N=NI-1 are put equal to zero,
    %     !     and the others are the same as described above. This means that an
    %     !     N'th order local covariance function will be computed. In this case,
    %     !     EPS must have NZ zero values stored.
    %     
    %     ! In all cases, N1 must be less than 300 and EPS must have dimension N1.

    % Constants
    RE = 6371.0; % Earth radius in km
    GM = 3.986004418e14; % Gravitational constant * Earth mass (m^3/s^2)
    RBJ2 = RE * 2;
    RBJ = sqrt(RBJ2);
    AM = RE / 5;
    AM2 = AM / 5;
 %  A IS IN UNITS OF MGAL**2, AM IN UNITS OF MGAL*M/SEC AND AM2 IN 
 %  UNITS OF (M/SEC)**2. RBJ IS THE RADIUS OF THE BJERHAMMAK-SPHERE. 
    MODEL = (N1 == 1);
    EPSC = zeros(1, 300);
 %  WE WILL NOW COMPUTE THE MODIFIED (POTENTIAL) DEGREE-VARIANCES, CF. 
 %  EQUATION (151)
    if ~MODEL && N1 >= 3
        for I = 3:N1
            RI = double(I);
            if I == 3
                EPS(3) = EPS(3) * RBJ2 * 1.0e-10;
            else
                EPS(I) = RBJ2 * (EPS(I) / ((RI - 1) ^ 2)) * 1.0e-10 - AM / ((RI - 1) * (RI - 2) * (RI + 8));
            end
        end
    end
end

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
        EPSC(I) = EPS(I) * (((I - 2) * 5 / RBJ) ^ IP);
    end
    NOTD = (KTYPE == 1) || (KTYPE == 3.3) || (KTYPE == 3.7);
    NOTDD = (KTYPE ~= 5) && (KTYPE ~= 4);
end

function COVC(PSI,HP,HQ,COV) 
    %  BY THE CALL OF COVC THE COVARIANCE OF TYPE KTYPE WILL RE COMPUTED FOK 
    %  PCINTS P AND 0 HAVING SPHERICAL DISTANCE (RADIANS) PSI, WHERE HP IS 
    %  THE HEIGHT OF P AHiJVE THE EARTH AND HQ THE HEIGHT OF Q ABOVE THE 
    %  EARTH, THE CBVARIANCE WILL RE RETURNED BY THE VARIABLE CBV, UNITS ARE 
    %  PRODUCTS OF MGAL, METERS AND ARCSECONDS.
    T = cos(PSI);
    U = sin(PSI);
    T2 = T * T;
    U2 = U * U;
    RP = RE + HP;
    RQ = RE + HQ;
    S = RBJ2 / (RP * RQ);
    S2 = S * S;
    S3 = S2 * S;
    TS = T * S;
    P2 = (D3 * D2 - D1) / D2;
    GP = GM / (RP * RP);
    GQ = GM / (RQ * RQ);

    %  THF QUANTITIES L,M AND N DFFINED IN EQ.(75) ARE HERE CALLFI) SL,SM 
    %  AND SN. L**2 = SL2. 
    SL2 = nl+S7-D2*TS 
    SL = DSORT(SL2i 
    SL3 = SL2*SL 
    SN = 01-TS+SL 
    SM = D1-TS-SL 
    SLN = SL*SN 
    SLNL = -DLOG( SN/D2 
    % WHEN WE ARE COMPUTING A LOCAL N'TH ORDER COVARIANCF OR A COVARIANCE 
    % FROM A GLORAL MODEL WITH EMPIRICAL DEGREE-VARIANCES UP TO AND INCLUSIVE DEGREE N, WE WILL HAVE TO COMPUTE THE SUM (154), THE SUM (155) 
    % (WHEN NOTD IS FALSE) AND THE SUM (156) (WHEN NOTDD IS FALSE). (154) 
    % WILL BE ACCUMMULATEO IN BO, (155) IN DBO AND (156) IN DDBO. 
    % WHEN THE VARIABLE MODEL IS TRUE, BO, DBO AND DDB0 WILL BE PUT EQUAL 
    % TO ZERO.
       
        BO = DO;
        DBO = DO;
        DDBO = DO;
        
        if ~MODEL
            B0 = 0;
            B1 = 0;
            B2 = 0;
            DB1 = 0;
            DDB1 = 0;
            RL1 = double(N1);
    % WE WILL NOW USE THE RECURSION FORMULAE (183),(185) AND (186) WHERE 
    % THE TERM (176A) DIVIDED BY T IS CALLED EL AND FL1 IS THE TERM (1768) 
    % FOR SUBSCRIPT L+l. 
            for I = 1:N1
                EL = (2 * RL1 - 1) * S / RL1;
                FL1 = -RL1 * S2 / (RL1 + 1);
                RL1 = RL1 - 1;
                
                B2 = B1;
                B1 = B0;
                B0 = B1 * EL * T + B2 * FL1 + EPSC(N1);
                
                DB2 = DB1;
                DB1 = DBO;
                DBO = EL * (DB1 * T + B1) + FL1 * DB2;
                
                DDB2 = DDB1;
                DDB1 = DDBO;
                DDBO = EL * (DB1 * 2 + DDB1 * T) + FL1 * DDB2;
            end
        end
        
        COV = B0 * S;
    end
    
    % Computation of closed expressions
    
    DPL = D1 + SL;
    DML = D1 - SL;
    P31 = D1 + TS + D1;
    B0  = B0 * S;
    FM1 = S * (SM + TS * SLNL);
    FM2 = S * (SM * P31 / D2 + S2 * (P2 * SLNL + U2 / D4));
    F1  = log(G1 + D2 * S / (D1 - S + SL));
    F2  = (SL - D1 + T * F1) / S;
    
    if ~NOTD
        DBO = DBO * S;
        
        DFM1 = S2 * (DML / SL + SLNL + TS * (D1 / SLN + D1 / SN));
        DFM2 = S2 * ((P31 / SL + D2 - 7.0 * TS - D3 * SL) / D2 + ...
             S * (D3 * T * SLNL + S * P2 * DPL / SLN));
        DF1  = S2 / SLN;
        DF2  = -D1 / SL + TS / SLN + F1 / S;
        DL   = -S / SL;
        
        if ~NOTDD
            DDBO = DDBO * S;
            
            DDFM1 = S3 * (D1 / SL3 + D2 * DPL / SLN + TS * ...
                   (D1 / (SL3 * SN) + (DPL / SLN) ^ 2));
            DDFM2 = S3 * ((6.0 / SL + P31 / SL3 - 7.0) / D2 + ...
                   D3 * SLNL + 6.0 * TS * (DPL / SLN) + ...
                   P2 * S2 * ((DPL / SLN) ^ 2 + D1 / (SL3 * SN)));
            DDF1  = S3 * (DPL / SLN ^ 2 + D1 / (SN * SL3));
            DDF2  = (-S2 / SL3 + D2 * DF1 + T * DDF1) / S;
            DDL   = -S2 / SL3;
        end
    end
    
    for I = 2:N1
        RI  = double(I);
        D12 = D2 * RI - D1;
        D11 = (RI - D1) / S;
        
        FB  = (SL + D12 * T * F2 - D11 * F1) / (RI * S);
        F1  = F2;
        F2  = FB;
        
        if ~NOTD
            DFB  = (DL + D12 * (F1 + T * DF2) - D11 * DF1) / (RI * S);
            DF1  = DF2;
            DF2  = DFB;
            
            if ~NOTDD
                DDFB  = (DDL + D12 * (D2 * DF1 + T * DDF2) - D11 * DDF1) / (RI * S);
                DDF1  = DDF2;
                DDF2  = DDFB;
            end
        end
    end
    
    if ~(NOTD || KTYPE == 2)
        DK = DBO + AM2 * RBJ2 * (IB1 * DFM2 - IB2 * DFM1 - D3 * T * S3) + ...
             DFB - S2 / IB1 - D3 * S3 * (1.0 / IB2) / IR12;
    end
    
    switch KTYPE
        case 1
            COV = S * RO + A * S * (IR1 * (FB - S / B - S2 * T / IR2) + FM2) / IB2;
        case 2
            COV = D * (DBO * KRJ / (KP * KQ) + AM * S * K * ...
                 (DFM2 - DFB + S2 / IR1 + D3 * S3 * T / IR2) / IB2) / (G * S * KADSEC);
        case 3
            COV = (HO * RHJ + AM * RHJ * (FM2 - FB + S / B + ...
                 S2 * T / IR2 + S3 * P2 / IB2) / IR2) / (RP * GQ);
        case 4
            COV = (T - DK / (KP * KQ) - U2 * ...
                 (DDBO / (KP * KQ) + AM2 * S * (IR1 * DDF1 + ...
                 U2 - IR2 * (DDFM1 - D3 * S3) + DDFB - D3 * S3 / IR2) / IB12)) * ...
                 KADSEC / (GP * GQ);
        case 5
            COV = DK / (RP * RQ * GP * GQ) * RADSE2;
        case 6
            COV = U * DK / (GP * GQ * RP) * RADSEC;
        case 7
            COV = (DBO + AM2 * RBB2 * ...
                 (IB1 * FM2 - IB2 * (FM1 - S3 * P2) + ...
                 FB - S / R - S2 * T / IB1 - S3 * P2 / IB2) / IB12) / (GP * GQ);
end
