SUBROUTINE COVA(EPSrN1) 
C THE SllBROUTINE COMPUTES ONE OF SEVEN DIFFERENT COVARIANCES (SEE RELOW)? USING THE ANOMALY DEGREE-VARIANCE MODEL GIVEN THROUGH THt VAL- 
C (JES OF TABLE SEVEN AND EQUATION (68)e(THE QUANTITY S IN THE TABLE IS 
C HERE CALLED SE). 
C THERE ARE THREE ENTRIES TO THE SUBROUTINE, WHICH HAVE TO BE CALLED IN 
C THE SEQUENCE COVA?COVB AND COVC. 
C BY THE CALL OF COVA, THE KIND OF COVARIANCE FUNCTION TO BE USED IS 
C DETERMINED. THERE ARE THREE POSSIBILITIES: 
C (1) THE COVARIANCE MODEL FOUR (EQUATIONS (130)-(132) AND (136)-(139)) 
C IS USED WITHOUT MODIFICATIONS. IN THIS CASE EPS WILL BE A DUMMY 
C ARRAY AND N1 MUST BE EQUAL TO ONE. 
C THE LOGICAL VARIABLE MODEL WILL GET THE VALUE TRUE IN THIS CASE, 
C (2) A NUMBER (NI) OF THE ANOMALY DEGREE-VARIANCES (DEGREE ZERO TO 
C NI-l) ARE PUT EQUAL TO EMPIRICAL DETERMINED DEGREE-VAKIANCES. 
C THE DEGREE-VARIANCE OF DEGRFE K WILL HAVE TO BE STOHFD IN 
C EPS(K+l) (IN UNITS OF MGAL**2). 
C 13) THE DEGREE-VARIANCES CIF DEGREE zEun TU N = NI-i ARE PUT EQUAL TO 
C ZEKOe(AND THE OTHERS ARE THE SAME AS ABOVE DESCK1HED)e THIS MEANS 
C THAT AN N'TH ORDER LOCAL COVARIANCE FUNCTION WILL BE COMPUTED. IN 
C THIS CASE EPS MUST HAVE NZ ZERO VALUES STORED, 
C IN ALL CASES N1 MUST BE LESS THAN 300 AND EPS MUST HAVE DIMENSION 
C N1. 
IMPLICIT REAL *8(A-HTO-Z) 
LOGICAL MODEL, NOTDV NOTDD 
DIMENSION EPSC(300)vEPS(l) 
DATA RE?GMyApSEIB, 1019 I829 IBMlyEpSC( ~)TEPSC(~)?DOT~~TD~TD~TD~~ 
*D5,RADSEC/6371~0D3~3~98014~425~28D0~0~999617D0~24~0D0~25~26~23~ 
*3*O. ODO,laODOt 2.0D09 3e0D09400D0v lmOD59 206264.806DO/ 
I012 = IB13kIB2 
KADSE2 = RADSEC**2 
RE2 = RE*RE 
RBJ2 = REZISE 
RBJ = DSQRT(RHJ2) 
AM = A/D5 
AM2 = AM/D5 
C A IS IN UNITS OF MGAL**2r AM IN UNITS OF MGAL*M/SEC AND AM2 IN 
C UNITS OF (M/SEC)zc*Z. RBJ IS THE RADIUS OF THE BJERHAMMAK-SPHERE. 
MODEL = N1 . EQ. 1
   IF ( MODEL) GO TO 20

   IF (N1 .LT. 3) GO TO 20
   DO 10 I = 3, N1
   RI = DFLOAT(I-1)
   IF (I .EQ. 3)EPS(3) = EPS(3) * RBJ2 * 1.0D-10
10 IF (I .GT. 3)EPS(I) = 
  *   RBJ2 * (EPS(I) / ((RI - D1) ** 2)) * 1.0D-10 - AM2 / ((RI - D1) * (RI - D2) * (RI + 8)))
 
20 RETURN 

 ENTRY COVB( KTYPE )
C BY THE CALL OF CClVB9 THE TYPE [IF COVARIANCE TO HF C0I"IPIJ'TkD IS I)FTFRMHNEU HY THF VALUE OF KTYPE, SO THAT WE GE1 THk CIlVAKIANCt HET 
C THE GRAVITY ANOMALY AT P AND 1HE GRAVITY ANOMALY AT (.l FOK KTYPE=1, 
C THE - - --- THE LONGITUOI ONAL COMPONtNT 
C OF THE DEFLECTION OF THE VFKTJCAL AT O FOR KTYtJE=2, 
C THE - - AT P AND THE HEIGHT ANOMALY AT O FOR KTYPEz3, 
C THE LONGITUDIONAL COMPONENT OF THE DEFLECTION OF THE VERTICAL AT P AND THE SAME TYPE OF QUANTITY AT Q FOR KTYPE=4, 
C THE TRANSVERSAL CCIMPOMENT OF THE DEFLECTION OF THE VFKTICAL AT P AND THE SAME TYPE OF QUANTITY AT fJ FOR KTYPEz59 
C THE LONGITIJOIONAL COMPONENT OF THE DEFLECTION OF THE VEKTICAL AT P AND THE HEIGHT ANOMALY AT Q FOR KTYPEz6, 
C AND THE HEIGHT ANOMALY AT P AND THE HEIGHT ANOMALY AT Q FOR KTYPE=7, 
C THE VkLUE OF KTYPE WILL THEN ALSO 9ETEKMINE WHICH OF THE COEFFICIENTS 
C (151)-(153)9THAT WE WILL USE IN THE EVALUATION OF THF LEGENDRE-SERIES 
C ANCI Whether NO [lIFFfRfNTIATI~lNT D1 FFEKENT IATI(1 ONE TIME OK DIFFEKENTIQTION TWO TIMES b!ITH RESPECT TCI THE VARIABLE T TAKFS PLACE. TWO 
C LOGICAL VAYIAHLES NOT11 AhlO NOTDI) AKE USED TO UISTINGIJISH BETWEEN THE 
C SITUATIONS. 
IF (MODEL) GO TO 35
C 
    IF (KTYPE.EQ.1) IP=2
    IF (KTYPE.EQ.2.OR.KTYPE.EQ.3) IP=1
    IF (KTYPE.GT.3) IP=0
    DO 30 I=3,N1
30 EPSC(I)=EPS(I)*((I-2)*D5/RBJ)**IP
35 NOTD = KTYPE.EQ. l,OR.KTYPE.E(3.3.OK.KTYPE,E(3.7 
   NOTDD = KTYPE,NE, 5, AND.KTYPE.NEe4 
   RETlJRN 
ENTRY COVC(PSI~WPPHO~COV) 
C BY THE CALL OF COVC THE COVARIANCE OF TYPE KTYPE WILL RE COMPUTED FOK 
C PCINTS P AND 0 HAVING SPHERICAL DISTANCE (RADIANS) PSI, WHERE HP IS 
C THE HEIGHT OF P AHiJVE THE EARTH AND HQ THE HEIGHT OF Q ABOVE THE 
C EARTH, THE CBVARIANCE WILL RE RETURNED BY THE VARIABLE CBV, UNITS ARE 
C PR(1DUCTS OF MGALy METERS AND AKCSECONDS.
T =DCOS(PSI)
U =DSIN(PSI)
T2 =T*T
U2 =U*U
RP =RE+HP
RQ =RE+HQ
S = RBJ2/(RP*RQ)
S2 = S*S
S3 = S2*S
TS = T*S
P2 =(D3*D2-D1)/D2
GP=GM/(RP*RP)
GQ=GM/(RQ*RQ)

C THF QUANTITIES LqM AND N DFFINED IN fQb(75) ARE HERE CALLFI) SLqSM 
C AND SN. L**2 = SL2, 
SL2 = nl+S7-D2*TS 
SL = DSORT(SL2i 
SL3 = SL2*SL 
SN = 01-TS+SL 
SM = D1-TS-SL 
SLN = SL*SN 
SLNL = -DLOG( SN/D2 
C WHEN WE ARE COMPUTING A LOCAL NITH ORDER COVARIANCF OK A COVAKIANCE 
C FROM A GLORAL MODEL WITH EMPIRICAL DEGREE-VARIANCES UP TO AND INCLUSIVE DEGREE Nt WE WILL HAVE TO COMPUTE THE SlJM (1541, THE SUM (155) 
C (WHEN NOTO IS FALSE) AND THE SUK (156) (WHEN iQUTDD IS FALSE). (154) 
C WILL BE ACCUMMULATEO IN BOY (155) IN DBO AND (156) IN DDBOe 
C WHEN THE VARIABLE MODEL IS TRUE, 807 DRO AND DD00 WILL BE POT EQUAL 
C TO ZERO, 
BO = Do 
DBO = DO 
DDBO = DO 
IF (MODEL ) GO TO 45 

B1= DO
DB1= DO
DDB1= DO
L1= N1
RL1=DFLOAT(L1)

C WE WILL NOW USE THE RECURSION FORMULAE (183)9(185) AND (18619 WHERE 
C THE TERM (176A) DIVIDED BY T IS CALLED EL AN[) FLP IS THE TERM '(1768) 
C FOR SUBSCRIPT L+l, 
DO 40 I = 19 N1 
EL = (DZ*RLI-DP)*S/RLS 
FL1 = -RLl*S2/(RLl+D1) 
Rh1 = RL1-D1 
B2 = B1 
B1 = R0 
B0 = BI*EL*T+B2*FLl+EPSC(Ll) 
IF (NOTD) GO TO 40 
DB2 = DB1 
DBI = DRO 
DBO = EL*(DBl*T+B1)+FLl*DB2 
IF (NOTDD) GO TO 40 
D0B2 = DDB1 
DDBl = DDBO 
DDBO = EL* ( DBl*D2+DDB1*T )+FL1*DDB2 
40 L1 = L1-1 

C COMPUTATION [IF CLOSED EXPRESSIONS. FIKST SOME AUXI LLIAKY OIJANTIT IES. 
C FM1 IS THF QUANTITY (86)~ FM2 IS (8717 F1 IS (99) AND F2 IS (100) 
45 DPL = Dl+SL 
DML = (11-SL 
P31 = DWTS+Dl 
B0 = BO*S 
FM1 = S*( SM+TS*SLNL) 
FM2 = S*(SM*P31/D2+S2*(P2*SLNL+U2/D4)) 
F1 = DLOG(Gl+D2*S/(Dl-S+SL)) 
F2 = (SL-Dl+T*Fl )/S 
IF (NOTD) GO TO 48 
DBO = DBO*S 
C DFMl IS THE OIJANTITY (901, DFM2 IS (9.21, DF1 IS (101) AND LaF2 IS 
C (103). 
DFMl = SZ*(DML/SL+SLNL+TS*(Dl/SLRI+Dl/SN)I 
DFM2 = S2* ( ( P31/SL+D2-7.0DOltTS-D3*SL )/DZ+S* ( D3*T*SLNL * +S*P2*DPL/SLN) ) 
DF1 = S2/SLN 
DF2 = -DL/SL+TS/SLN+Fl/S 
DL =-S/SL 
IF (NOTDD)  GO TO 48 
DDBO = DDBO*S
C DDFM1 IS THE QUANTITY (91),DDFM2 is (93),DDF1 is (102)and DDF2 is (104).
DDFM1= S3*(D1/SL3+D2*DPL/SLN+TS*(D1/(SL3*SN)+(DPL/SLN)**2))
DDFM2 = S3*( ( b. ODO/SL+P3l/SL3-7. ODO )/D2+D3*SLNL+6.0DO*TS*CDPL/SLN * +P2*S2*( (DPL/SLN)**2+D1/( SL3*SN) )) 
DDF1 = S3*( @PL/SLN**2+Dl/( SN'gSL3f ) 
DDF2 = (-S2/SL3+02:*DFl+T*DDFl)/S 
DDL = -S2/SL3 
C WE CAN NOW {JSE THE RECURSION FORMULAE (96). (97) AND (98) FOR THE 
C CnMP[JTATIflN OF THE QUANTITY (731 CALLED FH AND ITS DERIVATIVES DFH 
C AND DDFB. 
48 DO 50 I = 2, IBMl 
RI = DFLOAT(1) 
D12 = D2:zRI-D1 
D11 = (RI-Dl)/S 
FB = (SL+DI2*T*F2-DI1*F1)/(RI*S) 
F1 = F2 
F2 = F0 
IF (NOTD) GO TO 50 
DFB = (DL+DI2*(Fl+T*DF2)-DIl*DFl)/(RI*S) 
DF1 = DF2 
DF2 = DFR 
IF (NOTDD) GO TO 50 
DDFR = (DDL+DI2*(D2*DF1+T*DDF2)-DI1~*DDF1)/(RI*S) 
DDFl = D0F2 
DDF2 = OOFB 
50 CONTINUE

IF (NOTD.OR.KTYPE.EO.2) GO TU 60 
C FROM (133) WE HAVE: 
DK = DBO+AM2*RBJ2*( IBl*DFM2-I52*~OFM1-O3*T*S3)+DFB-S2/IBl-~)3*S3*1~/ 
ht IB2)/IR12 
60 GO TO (61,62?63,64,65,66r67),KTYPE 
C EQUATION (132) AND (146) GIVES: 
61 COV = S*RO+A*S*(IR1*(FB-S/B-S2*T/IR2)+FM?)/IB2 
GO TO 70 
C FQUATION (139) AND (150) GIVES: 
62 COV = O*( DBO*KRJ/( KP*KQ)+AM*cS:K( DFM2-DFB+S2/ IR1+D3*:S3*cT/If32) / 182) / * G(S*KADSEC 
GO TO 70 
C EWQTION (131) AND (145) GIVFS: 
63 COV = (HO*RHJ+AM*RHJ~*(FM~-FF~+S/B+S~*T/IR~+S~*P~/~B~)/IR~)/ * ( RPI;:GQ ) 
GO TO 70 
C EOUATIGM (136) AND (lL7) GIVES: 
64 COV = (Tr-DK/( KP*KQ)-U2hk( DDRO/ (KP;:KQ i+AlYZ*S:::i IR1*D9Flu2-IR2*( DoFM1 
 -D3*S3 )+DDF3-03*S3/IR2 )/IB12) )*KADSE2/( GP*:GQI 
GO TO 70 
C EOUATION (137) AND (148) GIVES: 
65 COV = DK/(RP*RQ*GP*GO)*RADSE2 
GO TO 70 
C EQUATION (138) AND 6 149) GIVES: 
66 COV = U*DK/( GP*GQ*RP )hYRADSEC 
GQ TO 70 
C AND EQUATION (371, (130) AN0 (144) GIVES: 
67 COV = ~BO+AM2*RBB2*~IB11*FM2-IB2*(FM1-S3*P2~+FB-S/R-S2*T/IB1-S3*P2 .C 
.,- /IB2)/IB12)/(GP*GO) 
70 RETURN 
END



