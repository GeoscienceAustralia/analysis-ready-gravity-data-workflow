% this is how to test COVA

haversineDistance=(0:pi/180: pi);

for i=1:1
cov0 =   COVA(1, num2str(i), haversineDistance,0, 0);
cov500 = COVA(1, num2str(i), haversineDistance,500*10^3, 500*10^3);

figure
plot(rad2deg(haversineDistance),cov0,'r-')
hold on 
plot(rad2deg(haversineDistance),cov500,'r*')
hold on 
plot(rad2deg(haversineDistance),0*haversineDistance,'b.')
title(['covariance model ',num2str(i)])
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
            EPSC(I) = EPS(I) * (((I - 2) * D5 / RBJ) ^ IP);
        end
        NOTD = (KTYPE == 1) || (KTYPE == 3) || (KTYPE == 7);
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

    SL2 = D1 + S2 - D2 * TS;
    SL = sqrt(L2);
    SL3 = L2 * SL;
    SN = D1 - TS + SL;
    SM = D1 - TS - SL;
    SLN = SL * SN;
    SLNL = -log(SN / D2);

    % WHEN WE ARE COMPUTING A LOCAL N'TH ORDER COVARIANCF OR A COVARIANCE 
    % FROM A GLORAL MODEL WITH EMPIRICAL DEGREE-VARIANCES UP TO AND INCLUSIVE DEGREE N, WE WILL HAVE TO COMPUTE THE SUM (154), THE SUM (155) 
    % (WHEN NOTD IS FALSE) AND THE SUM (156) (WHEN NOTDD IS FALSE). (154) 
    % WILL BE ACCUMMULATEO IN BO, (155) IN DBO AND (156) IN DDBO. 
    % WHEN THE VARIABLE MODEL IS TRUE, BO, DBO AND DDB0 WILL BE PUT EQUAL 
    % TO ZERO.
       
    B0 = D0;
    DB0 = D0;
    DDB0 = D0;
        
    if MODEL
        return;
    end
            
    B1 = D0;
    DB1 = D0;
    DDB1 = D0;
    L1= N1;
    RL1 = double(L1);
    % WE WILL NOW USE THE RECURSION FORMULAE (183),(185) AND (186) WHERE 
    % THE TERM (176A) DIVIDED BY T IS CALLED EL AND FL1 IS THE TERM (1768) 
    % FOR SUBSCRIPT L+l. 
    for I = 1:N1
        EL = (D2 * RL1 - D1) * S / RL1;
        FL1 = -RL1 * S2 / (RL1 + D1);
        RL1 = RL1 - 1;
        
        B2 = B1;
        B1 = B0;
        B0 = B1 * EL * T + B2 * FL1 + EPSC(L1);
        
        if ~NOTD
            continue;
        end
    
        DB2 = DB1;
        DB1 = DB0;
        DB0 = EL * (DB1 * T + B1) + FL1 * DB2;
    
        if ~NOTDD
            continue;
        end
    
        DDB2 = DDB1;
        DDB1 = DDB0;
        DDB0 = EL * (DB1 * D2 + DDB1 * T) + FL1 * DDB2;
    
        L1 = L1 - 1;
    end

    
    % COMPUTATION OF CLOSED EXPRESSIONS. FIRST SOME AUXILLIARY QUUANTITIES. 
    % FM1 IS THF QUANTITY (86), FM2 IS (87), F1 IS (99) AND F2 IS (100)
    
    DPL = D1 + SL;
    DML = D1 - SL;
    P31 = D3 + TS + D1;
    B0  = B0 * S;
    FM1 = S * (SM + TS * SLNL);
    FM2 = S * (SM * P31 / D2 + S2 * (P2 * SLNL + U2 / D4));
    F1  = log(D1 + D2 * S / (D1 - S + SL));
    F2  = (SL - D1 + T * F1) / S;
    if ~NOTD
        continue;
    end
    DB0 = DB0 * S;
    % DFMl IS THE OIJANTITY (90) DFM2 IS (92) DF1 IS (101) AND DF2 IS 
    % (103). 
        
    DFM1 = S2 * (DML / SL + SLNL + TS * (D1 / SLN + D1 / SN));
    DFM2 = S2 * ((P31 / SL + D2 - 7.0 * TS - D3 * SL) / D2 + ...
         S * (D3 * T * SLNL + S * P2 * DPL / SLN));
    DF1  = S2 / SLN;
    DF2  = -D1 / SL + TS / SLN + F1 / S;
    DL   = -S / SL;
    if ~NOTDD
        continue;
    end
    DDB0 = DDB0 * S;  
    %DDFM1 IS THE QUANTITY (91),DDFM2 is (93),DDF1 is (102)and DDF2 is (104).
    DDFM1 = S3 * (D1 / SL3 + D2 * DPL / SLN + TS * ...
           (D1 / (SL3 * SN) + (DPL / SLN) ^ 2));
    DDFM2 = S3 * ((6.0 / SL + P31 / SL3 - 7.0) / D2 + ...
           D3 * SLNL + 6.0 * TS * DPL / SLN + ...
           P2 * S2 * ((DPL / SLN) ^ 2 + D1 / (SL3 * SN)));
    DDF1  = S3 * (DPL / SLN ^ 2 + D1 / (SN * SL3));
    DDF2  = (-S2 / SL3 + D2 * DF1 + T * DDF1) / S;
    DDL   = -S2 / SL3;
    % C WE CAN NOW USE THE RECURSION FORMULAE (96), (97) AND (98) FOR THE 
    % C COMPUTATION OF THE QUANTITY (73) CALLED FB AND ITS DERIVATIVES DFB 
    % C AND DDFB. 
        
    for I = 2:IBM1
        RI  = double(I);
        DI2 = D2 * RI - D1;
        DI1 = (RI - D1) / S;
        
        FB  = (SL + DI2 * T * F2 - DI1 * F1) / (RI * S);
        F1  = F2;
        F2  = FB;
        
        if ~NOTD
            continue;
        end
    
        DFB  = (DL + DI2 * (F1 + T * DF2) - DI1 * DF1) / (RI * S);
        DF1  = DF2;
        DF2  = DFB;
            
        if ~NOTDD
            continue;
        end
    
        DDFB  = (DDL + DI2 * (D2 * DF1 + T * DDF2) - DI1 * DDF1) / (RI * S);
        DDF1  = DDF2;
        DDF2  = DDFB;
         
    end
    
        
    if ~NOTD || (KTYPE == 2)
        % From equation (133), we have:
        DK = DBO + AM2 * RBJ2 * (IB1 * DFM2 - IB2 * (DFM1 - D3 * T * S3) + ...
             DFB - S2 / IB1 - D3 * S3 * T / IB2) / IB12;
    end
        
    switch KTYPE
            % EQUATION (132) AND (146) GIVES: 
        case 1
            COV = S * BO + A * S * (IB1 * (FB - S / B - S2 * T / IB1- S3 * P2 / IB2) + FM2) / IB2;
            % EQUATION (139) AND (150) GIVES: 
        case 2
            COV = U * (DBO * RBJ / (RP * RQ) + AM * S * ...
                 (DFM2 - DFB + S2 / IB1 + D3 * S3 * T / IB2) / IB2) / (GQ * RADSEC);
        % EQUAION (131) AND (145) GIVES:
        case 3
            COV = (BO * RBJ + AM * RBJ2 * (FM2 - FB + S / B + ...
                 S2 * T / IB1 + S3 * P2 / IB2) / IB2) / (RP * GQ);
            % EOUATIGM (136) AND (l47) GIVES:
        case 4
            COV = (T * DK / (RP * RQ) - U2 * ...
                 (DDB0 / (RP * RQ) + AM2 * S * (IB1 * DDFM2 - ...
                  IB2 * (DDFM1 - D3 * S3) + DDFB - D3 * S3 / IB2) / IB12)) * ...
                 RADSE2 / (GP * GQ);
            % EOUATION (137) AND (148) GIVES: 
        case 5
            COV = DK / (RP * RQ * GP * GQ) * RADSE2;
            % EQUATION (138) AND (149) GIVES:
        case 6
            COV = U * DK / (GP * GQ * RP) * RADSEC;
            % AND EQUATION (37), (130) AND (144) GIVES:
        case 7
            COV = (B0 + AM2 * RBJ2 * ...
                 (IB1 * FM2 - IB2 * (FM1 - S3 * P2) + ...
                 FB - S / B - S2 * T / IB1 - S3 * P2 / IB2) / IB12) / (GP * GQ);
    end
end
