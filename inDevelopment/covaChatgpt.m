SUBROUTINE COVA(EPSrN1)
! This subroutine computes one of seven different covariances (see below)
! using the anomaly degree-variance model given through the values of 
! Table Seven and Equation (68). The quantity S in the table is 
! here called SE.

! There are three entries to the subroutine, which have to be called in
! the sequence COVA, COVB, and COVC.

! By calling COVA, the kind of covariance function to be used is determined.
! There are three possibilities:
! (1) The covariance model four (Equations (130)-(132) and (136)-(139))
!     is used without modifications. In this case, EPS will be a dummy
!     array and N1 must be equal to one. The logical variable MODEL will
!     get the value TRUE in this case.
! (2) A number (NI) of the anomaly degree-variances (degree zero to NI-1)
!     are put equal to empirically determined degree-variances. The degree-
!     variance of degree K will have to be stored in EPS(K+1) (in units
!     of MGAL**2).
! (3) The degree-variances of degree zero to N=NI-1 are put equal to zero,
!     and the others are the same as described above. This means that an
!     N'th order local covariance function will be computed. In this case,
!     EPS must have NZ zero values stored.

! In all cases, N1 must be less than 300 and EPS must have dimension N1.

IMPLICIT REAL*8 (A-H, O-Z)
LOGICAL MODEL, NOTDV, NOTDD
DIMENSION EPSC(300), EPS(1)

! Define constants and initial values
DATA RE, GMyApSEIB /1019, 829, IBMlyEpSC(1), TEPSC(1), DOT / ! Define initial parameters
RADSEC = 6371.0D0
AM = A / 5
AM2 = AM / 5
RE2 = RE * RE
RBJ2 = RE * 2
RBJ = DSQRT(RBJ2)

! Calculate radius of Bjerhammak-Sphere
MODEL = N1 - EQ. 1

! If model flag is false, skip to next section
IF (.NOT. MODEL) GO TO 20

! Compute the modified (potential) degree-variances using Equation (151)
IF (N1 .LT. 3) GO TO 20

DO 10 I = 3, N1
    RI = DFLOAT(I)
    IF (I .EQ. 3) THEN
        EPS(3) = EPS(3) * RBJ2 * 1.0D-10
    END IF

    IF (I .GT. 3) THEN
        EPS(I) = RBJ2 * (EPS(I) / ((RI - D1) ** 2)) * 1.0D-10 - A / ((RI - D1) * (RI - D2) * (RI + 8))
    END IF
10 CONTINUE

20 RETURN
END SUBROUTINE COVA

ENTRY COVB(KTYPE)
! Determines the type of covariance function to be computed based on KTYPE.
! The covariance function is set as follows:
! - KTYPE = 1: Covariance between the gravity anomaly at P and at Q.
! - KTYPE = 2: Covariance between the longitudinal component of the deflection 
!              of the vertical at P and Q.
! - KTYPE = 3: Covariance between the gravity anomaly at P and the height 
!              anomaly at Q.
! - KTYPE = 4: Covariance between the longitudinal component of the deflection 
!              of the vertical at P and the same at Q.
! - KTYPE = 5: Covariance between the transversal component of the deflection 
!              of the vertical at P and the same at Q.
! - KTYPE = 6: Covariance between the longitudinal component of the deflection 
!              of the vertical at P and the height anomaly at Q.
! - KTYPE = 7: Covariance between the height anomaly at P and at Q.
!
! The value of KTYPE also determines which coefficients (151)-(153) are used
! in the evaluation of the Legendre series and whether differentiation 
! (one or two times) with respect to the variable T is required.
! Logical variables NOTD and NOTDD are used to distinguish between these cases.

IF (MODEL) GO TO 35

! Set the exponent IP based on KTYPE
IF (KTYPE .EQ. 1) IP = 2
IF (KTYPE .EQ. 2 .OR. KTYPE .EQ. 3) IP = 1
IF (KTYPE .GT. 3) IP = 0

! Modify EPSC values based on IP
DO 30 I = 3, N1
    EPSC(I) = EPS(I) * ((I - 2) * D5 / RBJ) ** IP
30 CONTINUE

35 CONTINUE

! Define logical variables for differentiation conditions
NOTD = (KTYPE .EQ. 1) .OR. (KTYPE .EQ. 3) .OR. (KTYPE .EQ. 7)
NOTDD = (KTYPE .NE. 5) .AND. (KTYPE .NE. 4)

RETURN

ENTRY COVC(PSI, HP, HQ, COV)
! Computes the covariance of type KTYPE for points P and Q
! PSI  - Spherical distance (in radians) between P and Q
! HP   - Height of P above the Earth
! HQ   - Height of Q above the Earth
! COV  - Computed covariance (returned value)
!
! The covariance is returned in units of products of MGAL, meters, and arcseconds.

T  = DCOS(PSI)
U  = DSIN(PSI)
T2 = T * T
U2 = U * U

RP = RE + HP
RQ = RE + HQ

S  = RBJ2 / (RP * RQ)
S2 = S * S
S3 = S2 * S
TS = T * S

P2 = (D3 * D2 - D1) / D2
GP = GM / (RP * RP)
GQ = GM / (RQ * RQ)

RETURN

! Compute the quantities L*M and N as defined in Eq. (75), referred to here as SL*SM and SN.
! L**2 is represented as SL2.
SL2 = D1 + S2 - D2 * TS
SL  = DSQRT(SL2)
SL3 = SL2 * SL

SN  = D1 - TS + SL
SM  = D1 - TS - SL
SLN = SL * SN
SLNL = -DLOG(SN / D2)

! When computing a local N-th order covariance function or using a global model
! with empirical degree-variances up to and including degree N1, we need to compute:
!  - Sum (154), stored in BO
!  - Sum (155), stored in DBO (if NOTD is false)
!  - Sum (156), stored in DDBO (if NOTDD is false)
!
! When the variable MODEL is true, BO, DBO, and DDBO will be set to zero.

BO   = 0.0D0
DBO  = 0.0D0
DDBO = 0.0D0

IF (MODEL) GO TO 45

B1   = 0.0D0
DB1  = 0.0D0
DDB1 = 0.0D0

L1  = N1
RL1 = DFLOAT(L1)

! Using the recursion formulas (183), (185), and (186), 
! where EL corresponds to term (176A) divided by T,
! and FL1 is term (176B) for subscript L+1.

DO 40 I = 1, N1
    EL  = (D2 * RL1 - D1) * S / RL1
    FL1 = -RL1 * S2 / (RL1 + D1)
    RL1 = RL1 - D1

    B2 = B1
    B1 = B0
    B0 = B1 * EL * T + B2 * FL1 + EPSC(L1)

    IF (NOTD) GO TO 40

    DB2 = DB1
    DB1 = DBO
    DBO = EL * (DB1 * T + B1) + FL1 * DB2

    IF (NOTDD) GO TO 40

    DDB2 = DDB1
    DDB1 = DDBO
    DDBO = EL * (DB1 * D2 + DDB1 * T) + FL1 * DDB2

40  L1 = L1 - 1

45 CONTINUE
RETURN

! Computation of closed expressions. First, define auxiliary quantities.
! FM1 corresponds to Eq. (86), FM2 to Eq. (87),
! F1 to Eq. (99), and F2 to Eq. (100).

45  DPL = D1 + SL
    DML = D1 - SL
    P31 = D1 + TS + D1
    B0  = B0 * S
    FM1 = S * (SM + TS * SLNL)
    FM2 = S * (SM * P31 / D2 + S2 * (P2 * SLNL + U2 / D4))
    F1  = DLOG(G1 + D2 * S / (D1 - S + SL))
    F2  = (SL - D1 + T * F1) / S

    IF (NOTD) GO TO 48

    DBO = DBO * S

    ! DFM1 corresponds to Eq. (90), DFM2 to Eq. (92),
    ! DF1 to Eq. (101), and DF2 to Eq. (103).

    DFM1 = S2 * (DML / SL + SLNL + TS * (D1 / SLN + D1 / SN))
    DFM2 = S2 * ((P31 / SL + D2 - 7.0D0 * TS - D3 * SL) / D2 + 
         S * (D3 * T * SLNL + S * P2 * DPL / SLN))
    DF1  = S2 / SLN
    DF2  = -D1 / SL + TS / SLN + F1 / S
    DL   = -S / SL

    IF (NOTDD) GO TO 48

    DDBO = DDBO * S

    ! DDFM1 corresponds to Eq. (91), DDFM2 to Eq. (93),
    ! DDF1 to Eq. (102), and DDF2 to Eq. (104).

    DDFM1 = S3 * (D1 / SL3 + D2 * DPL / SLN + TS * 
           (D1 / (SL3 * SN) + (DPL / SLN) ** 2))
    DDFM2 = S3 * ((6.0D0 / SL + P31 / SL3 - 7.0D0) / D2 + 
           D3 * SLNL + 6.0D0 * TS * (DPL / SLN) + 
           P2 * S2 * ((DPL / SLN) ** 2 + D1 / (SL3 * SN)))
    DDF1  = S3 * (DPL / SLN ** 2 + D1 / (SN * SL3))
    DDF2  = (-S2 / SL3 + D2 * DF1 + T * DDF1) / S
    DDL   = -S2 / SL3

    ! Using recursion formulas (96), (97), and (98) 
    ! for computing FH and its derivatives DFH and DDFH.

48  DO 50 I = 2, N1
        RI  = DFLOAT(I)
        D12 = D2 * RI - D1
        D11 = (RI - D1) / S

        FB  = (SL + D12 * T * F2 - D11 * F1) / (RI * S)
        F1  = F2
        F2  = FB

        IF (NOTD) GO TO 50

        DFB  = (DL + D12 * (F1 + T * DF2) - D11 * DF1) / (RI * S)
        DF1  = DF2
        DF2  = DFB

        IF (NOTDD) GO TO 50

        DDFB  = (DDL + D12 * (D2 * DF1 + T * DDF2) - D11 * DDF1) / (RI * S)
        DDF1  = DDF2
        DDF2  = DDFB

50  CONTINUE

IF (NOTD .OR. KTYPE .EQ. 2) GO TO 60 

! From Equation (133), we have:
DK = DBO + AM2 * RBJ2 * (IB1 * DFM2 - IB2 * DFM1 - D3 * T * S3) + 
     DFB - S2 / IB1 - D3 * S3 * (1.0 / IB2) / IR12 

60  GO TO (61, 62, 63, 64, 65, 66, 67), KTYPE 

! Equation (132) and (146) give:
61  COV = S * RO + A * S * (IR1 * (FB - S / B - S2 * T / IR2) + FM2) / IB2 
    GO TO 70 

! Equation (139) and (150) give:
62  COV = D * (DBO * KRJ / (KP * KQ) + AM * S * K * 
     (DFM2 - DFB + S2 / IR1 + D3 * S3 * T / IR2) / IB2) / (G * S * KADSEC) 
    GO TO 70 

! Equation (131) and (145) give:
63  COV = (HO * RHJ + AM * RHJ * (FM2 - FB + S / B + 
     S2 * T / IR2 + S3 * P2 / IB2) / IR2) / (RP * GQ) 
    GO TO 70 

! Equation (136) and (147) give:
64  COV = (T - DK / (KP * KQ) - U2 * 
     (DDBO / (KP * KQ) + AM2 * S * (IR1 * DDF1 + 
     U2 - IR2 * (DDFM1 - D3 * S3) + DDFB - D3 * S3 / IR2) / IB12)) * 
     KADSEC / (GP * GQ) 
    GO TO 70 

! Equation (137) and (148) give:
65  COV = DK / (RP * RQ * GP * GQ) * RADSE2 
    GO TO 70 

! Equation (138) and (149) give:
66  COV = U * DK / (GP * GQ * RP) * RADSEC 
    GO TO 70 

! Equations (37), (130), and (144) give:
67  COV = (DBO + AM2 * RBB2 * 
     (IB1 * FM2 - IB2 * (FM1 - S3 * P2) + 
     FB - S / R - S2 * T / IB1 - S3 * P2 / IB2) / IB12) / (GP * GQ) 

70  RETURN 
END






