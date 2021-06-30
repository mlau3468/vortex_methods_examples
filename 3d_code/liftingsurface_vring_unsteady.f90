SUBROUTINE VORTEX(X,Y,Z,X1,Y1,Z1,X2,Y2,Z2,GAMA,U,V,W)
    !SUBROUTINE VORTEX CALCULATES THE INDUCED VELOCITY (U,V,W) AT A POI
    !(X,Y,Z) DUE TO A VORTEX ELEMENT VITH STRENGTH GAMA PER UNIT LENGTH
    !POINTING TO THE DIRECTION (X2,Y2,Z2)-(X1,Y1,Z1).
    PAY=3.141592654
    RCUT=1.0E-10
    !CALCULATION OF R1 X R2
    R1R2X=(Y-Y1)*(Z-Z2)-(Z-Z1)*(Y-Y2)
    R1R2Y=-((X-X1)*(Z-Z2)-(Z-Z1)*(X-X2))
    R1R2Z=(X-X1)*(Y-Y2)-(Y-Y1)*(X-X2)
    !CALCULATION OF (R1 X R2 )**2
    SQUARE=R1R2X*R1R2X+R1R2Y*R1R2Y+R1R2Z*R1R2Z
    !CALCULATION OF R0(R1/R(R1)-R2/R(R2))
    R1=SQRT((X-X1)*(X-X1)+(Y-Y1)*(Y-Y1)+(Z-Z1)*(Z-Z1))
    R2=SQRT((X-X2)*(X-X2)+(Y-Y2)*(Y-Y2)+(Z-Z2)*(Z-Z2))
    IF((R1.LT.RCUT).OR.(R2.LT.RCUT).OR.(SQUARE.LT.RCUT)) GOTO 1
    R0R1=(X2-X1)*(X-X1)+(Y2-Y1)*(Y-Y1)+(Z2-Z1)*(Z-Z1)
    R0R2=(X2-X1)*(X-X2)+(Y2-Y1)*(Y-Y2)+(Z2-Z1)*(Z-Z2)
    COEF=GAMA/(4.0*PAY*SQUARE)*(R0R1/R1-R0R2/R2)
    U=R1R2X*COEF
    V=R1R2Y*COEF
    W=R1R2Z*COEF
    GOTO 2
    !WHEN POINT (X,Y,Z) LIES ON VORTEX ELEMENT; ITS INDUCED VELOCITY IS
1   U=0.
    V=0.
    W=0.
2   CONTINUE
    RETURN
END


SUBROUTINE DECOMP(N,NDIM,A,IP)
    REAL A(NDIM,NDIM),T
    INTEGER IP(NDIM)
    !MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
    !N = ORDER OF MATRIX. NDIM = DECLARED DIMENSION OF ARRAY A.
    !A = MATRIX TO BE TRIANGULARIZED.
    !IP(K) , K .LT. N = INDEX OF K-TH PIVOT ROW.
    IP(N) = 1
    DO 6 K = 1, N
    IF(K.EQ.N) GOTO 5
    KP1 = K + 1
    M = K
    DO 1 I = KP1, N
    IF( ABS(A(I,K)).GT.ABS(A(M,K))) M=I
1   CONTINUE
    IP(K) = M
    IF(M.NE.K) IP(N) = -IP(N)
    T = A(M,K)
    A(M,K) = A(K,K)
    A(K,K) = T
    IF(T.EQ.0.E0) GO TO 5
    DO 2 I = KP1, N
2   A(I,K) = -A(I,K)/T
    DO 4 J = KP1, N
    T = A(M,J)
    A(M,J) = A(K,J)
    A(K,J) = T
    IF(T .EQ. 0.E0) GO TO 4
    DO 3 I = KP1, N
3   A(I,J) = A(I,J) + A(I,K)*T
4   CONTINUE
5   IF(A(K,K) .EQ. 0.E0) IP(N) = 0
6   CONTINUE
    RETURN
END

SUBROUTINE SOLVER(N,NDIM,A,B,IP)
    REAL A(NDIM,NDIM), B(NDIM), T
    INTEGER IP(NDIM)
    !SOLUTION OF LINEAR SYSTEM, A*X = B.
    !N = ORDER OF MATRIX.
    !NDIM = DECLARED DIMENSION OF THE ARRAY A.
    !B = RIGHT HAND SIDE VECTOR.
    !IP = PIVOT VECTOR OBTAINED FROM SUBROUTINE DECOMP.
    !B = SOLUTION VECTOR, X.
    IF(N.EQ.1) GOTO 9
    NM1 = N - 1
    DO 7 K = 1, NM1
    KP1 = K + 1
    M = IP(K)
    T = B(M)
    B(M) = B(K)
    B(K) = T
    DO 7 I = KP1, N
7    B(I) = B(I) + A(I,K)*T
    DO 8 KB = 1, NM1
    KM1 = N - KB
    K = KM1 + 1
    B(K) = B(K)/A(K,K)
    T = -B(K)
    DO 8 I = 1, KM1
8    B(I) = B(I) + A(I,K)*T
9    B(1) = B(1)/A(1,1)
    RETURN
END


SUBROUTINE WAKE(X,Y,Z,IT,U,V,W)
DIMENSION VORTIC(50,13),QW(50,14,3)
!CALCULATES SEMI WAKE INDUCED VELOCITY AT POINT (X,Y,Z) AT T=IT*DT,
!IN THE INERTIAL FRAME OF REFERENCE
COMMON VORTIC,QW
COMMON/NO2/ IB,JB,CH,SIGN
COMMON/NO3/ IW
U=0
V=0
W=0
I1=IT-1
DO 1 I=1,I1
DO 1 J=1,JB
VORTEK=VORTIC(I,J)
CALL VORTEX(X,Y,Z,QW(I,J,1),QW(I,J,2),QW(I,J,3),QW(I+1,J,1),QW(I+1,J,2),QW(I+1,J,3),VORTEK,U1,V1,W1)
CALL VORTEX(X,Y,Z,QW(I+1,J,1),QW(I+1,J,2),QW(I+1,J,3),QW(I+1,J+1,1),QW(I+1,J+1,2),QW(I+1,J+1,3),VORTEK,U2,V2,W2)
CALL VORTEX(X,Y,Z,QW(I+1,J+1,1),QW(I+1,J+1,2),QW(I+1,J+1,3),QW(I,J+1,1),QW(I,J+1,2),QW(I,J+1,3),VORTEK,U3,V3,W3)
CALL VORTEX(X,Y,Z,QW(I,J+1,1),QW(I,J+1,2),QW(I,J+1,3),QW(I,J,1),QW(I,J,2),QW(I,J,3),VORTEK,U4,V4,W4)
U=U+U1+U2+U3+U4
V=V+V1+V2+V3+V4
W=W+W1+W2+W3+W4
1 CONTINUE
RETURN
END

SUBROUTINE VELOCE(X,Y,Z,U,V,W,IT,JS1,JS2)
    DIMENSION GAMA(4,13)
    COMMON/NO1/ SX,SZ,CS1,SN1,GAMA
    COMMON/NO2/ IB,JB,CH,SIGN
    !SUBROUTINE VELOCE CALCULATES INDUCED VELOCITIES DUE TO THE WING
    !AND ITS WAKES IN A POINT (X,Y,Z) GIVEN IN THE INERTIAL FRAME OF
    !REFERENCE.
    X1=(X-SX)*CS1+(Z-SZ)*SN1
    Y1=Y
    Z1=-(X-SX)*SN1+(Z-SZ)*CS1
    CALL WAKE(X,Y,Z,IT,U1,V1,W1)
    CALL WAKE(X,-Y,Z,IT,U2,V2,W2)
    CALL WING(X1,Y1,Z1,GAMA,U3,V3,W3)
    CALL WING(X1,-Y1,Z1,GAMA,U4,V4,W4)
    U33=CS1*(U3+U4)-SN1*(W3+W4)
    W33=SN1*(U3+U4)+CS1*(W3+W4)
    !INFLUENCE OF MIRROR IMAGE
    IF(CH.GT.100.0) GOTO 1
    X2=(X-SX)*CS1+(-Z-SZ)*SN1
    Z2=-(X-SX)*SN1+(-Z-SZ)*CS1
    CALL WAKE(X,Y,-Z,IT,U5,V5,W5)
    CALL WAKE(X,-Y,-Z,IT,U6,V6,W6)
    CALL WING(X2,Y1,Z2,GAMA,U7,V7,W7)
    CALL WING(X2,-Y1,Z2,GAMA,U8,V8,W8)
    U77=CS1*(U7+U8)-SN1*(W7+W8)
    W77=SN1*(U7+U8)+CS1*(W7+W8)
    GOTO 2
1   CONTINUE
    U5=0.0
    U6=0.0
    U77=0.0
    V5=0.0
    V6=0.0
    V7=0.0
    V8=0.0
    W5=0.0
    W6=0.0
    W77=0.0
2   CONTINUE
    !VELOCITIES MEASURED IN INERTIAL FRAME
    U=U1+U2+U33+U5+U6+U77
    V=V1-V2+V3-V4+V5-V6+V7-V8
    W=W1+W2+W33-W5-W6-W77
    RETURN
END

SUBROUTINE WING(X,Y,Z,GAMA,U,V,W)
    DIMENSION GAMA(4,13),QF(5,16,3),A1(4,13),VORTIC(50,13),QW(50,14,3)
    DIMENSION ALF(5),SNO(5),CSO(5),VORT1(50,13),QW1(50,14,3)
    COMMON VORTIC,QW,VORT1,QW1,QF,A1
    COMMON IT,ALF,SNO,CSO
    COMMON/NO2/ IB,JB,CH,SIGN
    !CALCULATES SEMI WING INDUCED VELOCITY AT A POINT (X,Y,Z) DUE TO WI
    !VORTICITY DISTRIBUTION GAMA(I,J) IN A WING FIXED COORDINATE SYSTE
    U=0
    V=0
    W=0
    DO 7 I=1,IB
    DO 7 J=1,JB
    CALL VORTEX(X,Y,Z,QF(I,J,1),QF(I,J,2),QF(I,J,3),QF(I,J+1,1),QF(I,J+1,2),QF(I,J+1,3),GAMA(I,J),U1,V1,W1)
    CALL VORTEX(X,Y,Z,QF(I,J+1,1),QF(I,J+1,2),QF(I,J+1,3),QF(I+1,J+1,1),QF(I+1,J+1,2),QF(I+1,J+1,3),GAMA(I,J),U2,V2,W2)
    CALL VORTEX(X,Y,Z,QF(I+1,J+1,1),QF(I+1,J+1,2),QF(I+1,J+1,3),QF(I+1,J,1),QF(I+1,J,2),QF(I+1,J,3),GAMA(I,J),U3,V3,W3)
    CALL VORTEX(X,Y,Z,QF(I+1,J,1),QF(I+1,J,2),QF(I+1,J,3),QF(I,J,1),QF(I,J,2),QF(I,J,3),GAMA(I,J),U4,V4,W4)
    U0=U1+U2+U3+U4
    V0=V1+V2+V3+V4
    W0=W1+W2+W3+W4
    A1(I,J)=U0*SNO(I)+W0*CSO(I)
    IF(SIGN.GE.1.0) A1(I,J)=U0*SNO(I)-W0*CSO(I)
    U=U+U0
    V=V+V0
    W=W+W0
7   CONTINUE
    RETURN
END

SUBROUTINE WINGL(X,Y,Z,GAMA,U,V,W)
    DIMENSION GAMA(4,13),QF(5,16,3),A1(4,13),VORTIC(50,13),QW(50,14,3)
    DIMENSION ALF(5),SNO(5),CSO(5),VORT1(50,13),QW1(50,14,3)
    COMMON VORTIC,QW,VORT1,QW1,QF,A1
    COMMON IT,ALF,SNO,CSO
    COMMON/NO2/ IB,JB,CH,SIGN
    !CALCULATES INDUCED VELOCITY AT A POINT (X,Y,Z) DUE TO LONGITUDINAL
    !VORTICITY DISTRIBUTION GAMAX(I,J) ONLY(SEMI-SPAN), IN A WING FIXED
    !COORDINATE SYSTEM + (T.E. UNSTEADY VORTEX).
    !** SERVES FOR INDUCED DRAG CALCULATION ONLY **

    U=0.
    V=0.
    W=0.
    DO 7 I=1,IB
    DO 7 J=1,JB
    CALL VORTEX(X,Y,Z,QF(I,J+1,1),QF(I,J+1,2),QF(I,J+1,3),QF(I+1,J+1,1),QF(I+1,J+1,2),QF(I+1,J+1,3),GAMA(I,J),U2,V2,W2)
    CALL VORTEX(X,Y,Z,QF(I+1,J,1),QF(I+1,J,2),QF(I+1,J,3),QF(I,J,1),QF(I,J,2),QF(I,J,3),GAMA(I,J),U4,V4,W4)
    U=U+U2+U4
    V=V+V2+V4
    W=W+W2+W4
7   CONTINUE
    !ADD INFLUENCE OF LATEST UNSTEADY WAKE ELEMENT:
    I=IB
    DO 8 J=1,JB
    CALL VORTEX(X,Y,Z,QF(I+1,J+1,1),QF(I+1,J+1,2),QF(I+1,J+1,3),QF(I+1,J,1),QF(I+1,J,2),QF(I+1,J,3),GAMA(I,J),U3,V3,W3)
    U=U+U3
    V=V+V3
    W=W+W3
8   CONTINUE
    RETURN
END


SUBROUTINE GEO(B,C,S,AR,IB,JB,DX,DY,DGAP,ALFA)
    DIMENSION BB(13),ALF(5),SN(5),CS(5),SNO(5),CSO(5)
    DIMENSION QF(5,16,3),QC(4,13,3),DS(4,13)
    DIMENSION VORTIC(50,13),QW(50,14,3),A1(4,13),ALAMDA(2)
    DIMENSION VORT1(50,13),QW1(50,14,3)
    COMMON VORTIC,QW,VORT1,QW1,QF,A1
    COMMON IT,ALF,SNO,CSO,BB,QC,DS,ALAMDA,DXW
    PAY=3.141592654
    !IB:NO. OF CHORDWISE BOXES,JB:NO. OF SPANWISE BOXES
    IB1=IB+1
    JB1=JB+1
    DO 2 I=1,IB1
    SN(I)=SIN(ALF(I))
2   CS(I)=COS(ALF(I))
    CTG1=TAN(PAY/2.-ALAMDA(1))
    CTG2=TAN(PAY/2.-ALAMDA(2))
    CTIP=C+B*(CTG2-CTG1)
    S=B*(C+CTIP)/2.
    AR=2.*B*B/S
    !WING FIXED VORTICES LOCATION ( QF(I,J,(X,Y,Z))...)
    BJ=0.
    DO 3 J=1,JB1
    IF(J.GT.1) BJ=BJ+BB(J-1)
    Z1=0.
    DC1=BJ*CTG1
    DC2=BJ*CTG2
    DX1=(C+DC2-DC1)/IB
    !DC1=LEADING EDGE X,
    !DC2=TRAILING EDGE X
    DO 1 I=1,IB
    QF(I,J,1)=DC1+DX1*(I-0.75)
    QF(I,J,2)=BJ
    QF(I,J,3)=Z1-0.25*DX1*SN(I)
1   Z1=Z1-DX1*SN(I)
    !THE FOLLOWING LINES ARE DUE TO WAKE DISTANCE FROM TRAILING EDGE
    QF(IB1,J,1)=C+DC2+DXW
    QF(IB1,J,2)=QF(IB,J,2)
3   QF(IB1,J,3)=Z1-DXW*SN(IB)
    !WING COLLOCATION POINTS
    DO 4 J=1,JB
    Z1=0.
    BJ=QF(1,J,2)+BB(J)/2.
    DC1=BJ*CTG1
    DC2=BJ*CTG2
    DX1=(C+DC2-DC1)/IB
    DO 4 I=1,IB
    QC(I,J,1)=DC1+DX1*(I-0.25)
    QC(I,J,2)=BJ
    QC(I,J,3)=Z1-0.75*DX1*SN(I)
    Z1=Z1-DX1*SN(I)
4   DS(I,J)=DX1*BB(J)
    !ROTATION OF WING POINTS DUE TO ALFA
    SN1=SIN(-ALFA)
    CS1=COS(-ALFA)
    DO 6 I=1,IB1
    DO 6 J=1,JB1
    QF1=QF(I,J,1)
    QF(I,J,1)=QF1*CS1-QF(I,J,3)*SN1
    QF(I,J,3)=QF1*SN1+QF(I,J,3)*CS1
    IF((I.EQ.IB1).OR.(J.GE.JB1)) GOTO 6
    QC1=QC(I,J,1)
    QC(I,J,1)=QC1*CS1-QC(I,J,3)*SN1
    QC(I,J,3)=QC1*SN1+QC(I,J,3)*CS1
6   CONTINUE
    RETURN
    END

!------------------------------------------------------------------------
!PROGRAM No. 16: UNSTEADY RECTANGULAR LIFTING SURFACE (VLM)

!THIS IS A 3-D LINEAR CODE FOR RECTANGULAR PLANFORMS (WITH GROUND EFFECT)
!IN UNSTEADY MOTION USING THE VORTEX LATTICE METHOD (BY JOE KATZ, 1975).
DIMENSION ALF(5),SNO(5),CSO(5),ALAM(4),GAMA1J(5)
DIMENSION QF(5,16,3),QC(4,13,3),BB(13),DLY(13)
DIMENSION GAMA(4,13),DL(4,13),DP(4,13),DS(4,13),DLT(4,13),DD(4,13)
DIMENSION A1(4,13),QW(50,14,3),VORTIC(50,13),UVW(50,14,3)
DIMENSION QW1(50,14,3),VORT1(50,13),US(13)
DIMENSION A(52,52),GAMA1(52),WW(52),DW(52),IP(52),ALAMDA(2)
DIMENSION WTS(4,13),X15(50),Y15(50),Y16(50),Y17(50),Z15(50)
COMMON VORTIC,QW,VORT1,QW1,QF,A1
COMMON IT,ALF,SNO,CSO,BB,QC,DS,ALAMDA,DXW
COMMON/NO1/ SX,SZ,CS1,SN1,GAMA
COMMON/NO2/ IB,JB,CH,SIGN
COMMON/NO3/ IW


! This is a comment line; it is ignored by the compiler
!print *, 'Hello, World!'
!MODES OF OPERATION
!1. STEADY STATE
!: SET DT=DX/VT*IB*10, AND NSTEPS=5
!2. SUDDEN ACCELERATION : SET DT=DX/VT/4. AND NSTEPS= UP TO 50
!3. HEAVING OSCILLATIONS: BH=HEAVING AMPL. OM= FREQ.
!4. PITCH OSCILLATIONS : OMEGA=FREQUENCY, TETA=MOMENTARY ANGLE
!5. FOR COMPUTATIONAL ECONOMY THE PARAMETER IW MIGHT BE USED
!NOTE; INDUCED DRAG CALCULATION INCREASES COMPUTATION TIME AND
!CAN BE DISCONNECTED.
!INPUT DATA
IB=4
JB=13
NSTEPS=50
NW=10
IPROG1=1
DO 100 IPROG=1,IPROG1
!NW - IS THE NUMBER OF (TIMEWISE) DEFORMING WAKE ELEMENTS.
!IB = NUMBER OF CHORDWISE PANELS, JB = NO. OF SPANWISE PANELS
!NSTEPS = NO. OF TIME STEPS.
PAY=3.141592654
RO=1.
BH=0.0
OM=0.0
VT=50.0
    C=1.
    B=2 !semi span
    DX=C/IB
    DY=B/JB
    CH=1000000.*C
    !C=CHORD; C=ROOT CHORD, B=SEMI SPAN, DX,DY=PANEL DIMENSIONS,
    !CH=GROUND CLEARANCE, VT=FAR FIELD VELOCITY.
    ALFA1=5.0
    ALFAO=0.0
    ALFA=(ALFA1+ALFAO)*PAY/180.0
    DO 2 I=1,IB
2   ALF(I)=0.
    ALF(IB+1)=ALF(IB)
    !ALF(IB+1) IS REQUIRED ONLY FOR QF(I,J,K) CALCULATION IN GEO.
    ALAMDA(1)=90.*PAY/180.
    ALAMDA(2)=ALAMDA(1)
    !ALAMDA(I) ARE SWEP BACK ANGLES. (ALAMDA < 90, SWEEP BACKWARD).
    DT=DX/VT/16
    T=-DT
    !TIME IN SECONDS
    DXW=0.3*VT*DT
    DO 3 J=1,JB
3   BB(J)=DY
    !CONSTANTS
    K=0
    DO 1 I=1,IB
    DO 1 J=1,JB
    K=K+1
    WW(K)=0.0
    DLT(I,J)=0.0
    VORTIC(I,J)=0.0
    VORT1(I,J)=0.0
    !GAMA(I,J)=1. IS REQUIRED FOR INFLUENCE MATRIX CALCULATIONS.
1   GAMA(I,J)=1.
!   CALCULATION OF COLLOCATION POINTS.

    CALL GEO(B,C,S,AR,IB,JB,DX,DY,0.0,ALFA)
    !GEO CALCULATES WING COLLOCATION POINTS QC,AND VORTEX TIPS QF
    WRITE(6,101)
    ALAM1=ALAMDA(1)*180./PAY
    ALAM2=ALAMDA(2)*180./PAY
    WRITE(6,102) ALFA1,ALAM1,B,C,ALAM2,S,AR,IB,JB,CH
    DO 31 I=1,IB
    ALL=ALF(I)*180./PAY
31    WRITE(6,111) I,ALL
    DO 4 I=1,JB,2
4    WRITE(6,105) I,BB(I)
    IB1=IB+1
    JB1=JB+1
    !=============
    !PROGRAM START
    !=============
    DO 100 IT=1,NSTEPS
    T=T+DT
    !PATH INFORMATION
    SX=-VT*T
    DSX=-VT
    CH1=CH
    IF(CH.GT.100.0) CH1=0.0
    SZ=BH*SIN(OM*T)+CH1
    DSZ=BH*OM*COS(OM*T)
    !DSX=DSX/DT
    DSZ=DSZ/DT
    TETA=0.0
    OMEGA=0.0
    VT=-COS(TETA)*DSX-SIN(TETA)*DSZ
    SN1=SIN(TETA)
    CS1=COS(TETA)
    WT=SN1*DSX-CS1*DSZ
    DO 6 I=1,IB
    SNO(I)=SIN(ALFA+ALF(I))
6   CSO(I)=COS(ALFA+ALF(I))
    !===========================
    !VORTEX WAKE SHEDDING POINTS
    !===========================
    DO 7 J=1,JB1
    QW(IT,J,1)=QF(IB1,J,1)*CS1-QF(IB1,J,3)*SN1+SX
    QW(IT,J,2)=QF(IB1,J,2)
    QW(IT,J,3)=QF(IB1,J,1)*SN1+QF(IB1,J,3)*CS1+SZ
7   CONTINUE
    !========================
    !AERODYNAMIC CALCULATIONS
    !========================'
    
    !INFLUENCE COEFFICIENTS CALCULATION
    K=0
    DO 14 I=1,IB
    DO 14 J=1,JB
    SIGN=0.0
    K=K+1
    IF(IT.GT.1) GOTO 12
    !MATRIX COEFFICIENTS CALCULATION OCCURS ONLY ONCE FOR THE
    !TIME-FIXED-GEOMETRY WING.
    CALL WING(QC(I,J,1),QC(I,J,2),QC(I,J,3),GAMA,U,V,W)
    L=0
    DO 10 I1=1,IB
    DO 10 J1=1,JB
    L=L+1
    !A(K,L) - IS THE NORMAL VELOCITY COMPONENT DUE TO A UNIT VORTEX
    !LATTICE.
10  A(K,L)=A1(I1,J1)
    !ADD INFLUENCE OF WING OTHER HALF PART
    CALL WING(QC(I,J,1),-QC(I,J,2),QC(I,J,3),GAMA,U,V,W)
    L=0
    DO 11 I1=1,IB
    DO 11 J1=1,JB
    L=L+1
11  A(K,L)=A(K,L)+A1(I1,J1)
    IF(CH.GT.100.0) GOTO 12
    !ADD INFLUENCE OF MIRROR IMAGE.
    SIGN=10.0
    XX1=QC(I,J,1)*CS1-QC(I,J,3)*SN1+SX
    ZZ1=QC(I,J,1)*SN1+QC(I,J,3)*CS1+SZ
    XX2=(XX1-SX)*CS1+(-ZZ1-SZ)*SN1
    ZZ2=-(XX1-SX)*SN1+(-ZZ1-SZ)*CS1
    CALL WING(XX2,QC(I,J,2),ZZ2,GAMA,U,V,W)
    L=0
    DO 8 I1=1,IB
    DO 8 J1=1,JB
    L=L+1
8   A(K,L)=A(K,L)+A1(I1,J1)
    !ADD MIRROR IMAGE INFLUENCE OF WING'S OTHER HALF.
    CALL WING(XX2,-QC(I,J,2),ZZ2,GAMA,U,V,W)
    L=0
    DO 9 I1=1,IB
    DO 9 J1=1,JB
    L=L+1
9   A(K,L)=A(K,L)+A1(I1,J1)
    SIGN=0.0
12  CONTINUE
    IF(IT.EQ.1) GOTO 13
    !CALCULATE WAKE INFLUENCE
    XX1=QC(I,J,1)*CS1-QC(I,J,3)*SN1+SX
    ZZ1=QC(I,J,1)*SN1+QC(I,J,3)*CS1+SZ
    CALL WAKE(XX1,QC(I,J,2),ZZ1,IT,U,V,W)
    CALL WAKE(XX1,-QC(I,J,2),ZZ1,IT,U1,V1,W1)
    IF(CH.GT.100) GOTO 121
    CALL WAKE(XX1,QC(I,J,2),-ZZ1,IT,U2,V2,W2)
    CALL WAKE(XX1,-QC(I,J,2),-ZZ1,IT,U3,V3,W3)
    GOTO 122
121 U2=0.0
    U3=0.0
    V2=0.0
    V3=0.0
    W2=0.0
    W3=0.0
122 CONTINUE
    !WAKE INDUCED VELOCITY IS GIVEN IN INERTIAL FRAME
    U=U+U1+U2+U3
    W=W+W1-W2-W3
    U11=U*CS1+W*SN1
    W11=-U*SN1+W*CS1
    !WW(K) IS THE PREPENDICULAR COMPONENT OF WAKE INFLUENCE TO WING.
    WW(K)=U11*SNO(I)+W11*CSO(I)
13  CONTINUE
    !CALCULATE WING GEOMETRICAL DOWNWASH
    DW(K)=-VT*SNO(I)+QC(I,J,1)*OMEGA-WT
    !FOR GENERAL MOTION DW(K)=-VT*SIN(ALFA)+OMEGA*X
    WTS(I,J)=W11
    !W11 - IS POSITIVE SINCE THE LATEST UNSTEADY WAKE ELEMENT IS
    !INCLUDED IN SUBROUTINE WING
14  CONTINUE
    !SOLUTION OF THE PROBLEM: DW(I)=WW(I)+A(I,J)*GAMA(I)
    K1=IB*JB
    DO 15 K=1,K1
15  GAMA1(K)=DW(K)-WW(K)
    IF(IT.GT.1) GOTO 16
    !FOR NONVARIABLE WING GEOMETRY (WITH TIME), MATRIX INVERSION
    !IS DONE ONLY ONCE.
    CALL DECOMP(K1,52,A,IP)
16  CONTINUE
    CALL SOLVER(K1,52,A,GAMA1,IP)
    !HERE
    !*
    !THE SAME ARRAY SIZE IS REQUIRED,
    !AS SPECIFIED IN THE BEGINNING OF THE CODE
    !WING VORTEX LATTICE LISTING
    K=0
    DO 17 I=1,IB
    DO 17 J=1,JB
    K=K+1
17  GAMA(I,J)=GAMA1(K)
    !WAKE SHEDDING
    DO 171 J=1,JB
    !LATEST WAKE ELEMENTS LISTING
162 VORTIC(IT,J)=GAMA(IB,J)
    VORTIC(IT+1,J)=0.0
171 CONTINUE
    !===========================
    !WAKE ROLLUP CALCULATION
    !===========================
    IW=1
    IF(IT.EQ.1) GOTO 193
    IF(IT.GE.NW) IW=IT-NW+1
    !NW IS THE NUMBER OF (TIMEWISE) DEFORMING WAKE ELEMENTS.
    I1=IT-1
    JS1=0
    JS2=0
    DO 18 I=IW,I1
    DO 18 J=1,JB1
    CALL VELOCE(QW(I,J,1),QW(I,J,2),QW(I,J,3),U,V,W,IT,JS1,JS2)
    UVW(I,J,1)=U*DT
    UVW(I,J,2)=V*DT
    UVW(I,J,3)=W*DT
18  CONTINUE
    DO 19 I=IW,I1
    DO 19 J=1,JB1
    QW(I,J,1)=QW(I,J,1)+UVW(I,J,1)
    QW(I,J,2)=QW(I,J,2)+UVW(I,J,2)
    QW(I,J,3)=QW(I,J,3)+UVW(I,J,3)
19  CONTINUE
193 CONTINUE
    !==================
    !FORCES CALCULATION
    !==================
    FL=0.
    FD=0.
    FM=0.
    FG=0.0
    QUE=0.5*RO*VT*VT
    DO 20 J=1,JB
    SIGMA=0.
    SIGMA1=0.0
    DLY(J)=0.
    DO 20 I=1,IB
    IF(I.EQ.1) GAMAIJ=GAMA(I,J)
    IF(I.GT.1) GAMAIJ=GAMA(I,J)-GAMA(I-1,J)
    DXM=(QF(I,J,1)+QF(I,J+1,1))/2.
    !DXM IS VORTEX DISTANCE FROM LEADING EDGE
    SIGMA1=(0.5*GAMAIJ+SIGMA)*DX
    SIGMA=GAMA(I,J)
    DFDT=(SIGMA1-DLT(I,J))/DT
    !DFDT
    !IS THE VELOCITY POTENTIAL TIME DERIVATIVE
    DLT(I,J)=SIGMA1
    DL(I,J)=RO*(VT*GAMAIJ+DFDT)*BB(J)*CSO(I)
    !INDUCED DRAG CALCULATION
    CALL WINGL(QC(I,J,1),QC(I,J,2),QC(I,J,3),GAMA,U1,V1,W1)
    CALL WINGL(QC(I,J,1),-QC(I,J,2),QC(I,J,3),GAMA,U2,V2,W2)
    IF(CH.GT.100.0) GOTO 194
    XX1=QC(I,J,1)*CS1-QC(I,J,3)*SN1+SX
    ZZ1=QC(I,J,1)*SN1+QC(I,J,3)*CS1+SZ
    XX2=(XX1-SX)*CS1+(-ZZ1-SZ)*SN1
    ZZ2=-(XX1-SX)*SN1+(-ZZ1-SZ)*CS1
    CALL WINGL(XX2,QC(I,J,2),ZZ2,GAMA,U3,V3,W3)
    CALL WINGL(XX2,-QC(I,J,2),ZZ2,GAMA,U4,V4,W4)
    GOTO 195
194 W3=0.
    W4=0.
195 W8=W1+W2-W3-W4
    !ADD INFLUENCE OF MIRROR IMAGE (GROUND).
    CTS=-(WTS(I,J)+W8)/VT
    DD1=RO*BB(J)*DFDT*SNO(I)
    DD2=RO*BB(J)*VT*GAMAIJ*CTS
    DD(I,J)=DD1+DD2
    DP(I,J)=DL(I,J)/DS(I,J)/QUE
    DLY(J)=DLY(J)+DL(I,J)
    FL=FL+DL(I,J)
    FD=FD+DD(I,J)
    FM=FM+DL(I,J)*DXM
    FG=FG+GAMAIJ*BB(J)
20  CONTINUE
    CL=FL/(QUE*S)
    CD=FD/(QUE*S)
    CM=FM/(QUE*S*C)
    CLOO=2.*PAY*ALFA/(1.+2./AR)
    IF(ABS(CLOO).LT.1.E-20) CLOO=CL
    CLT=CL/CLOO
    CFG=FG/(0.5*VT*S)/CLOO
    !======
    !OUTPUT
    !======
    !PLACE PLOTTER OUTPUT HERE (e.g. T,SX,SZ,CL,CD,CM)
    !OTHER OUTPUT
    WRITE(6,106) T,SX,SZ,VT,TETA,OMEGA
    WRITE(6,104) CL,FL,CM,CD,CLT,CFG
    I2=5
    print *, IT
    IF(IT.NE.I2) GOTO 100
    WRITE(6,110)
    DO 21 J=1,JB
    DO 211 I=2,IB
211 GAMA1J(I)=GAMA(I,J)-GAMA(I-1,J)
    DLYJ=DLY(J)/BB(J)
21  WRITE(6,103) J,DLYJ,DP(1,J),DP(2,J),DP(3,J),DP(4,J),GAMA(1,J), GAMA1J(2),GAMA1J(3),GAMA1J(4)
    IF(IT.NE.I2) GOTO 100
    WRITE(6,107)
    DO 23 I=1,IT
    WRITE(6,109) I,(VORTIC(I,K1),K1=1,13)
    DO 23 J=1,3
    WRITE(6,108) J,(QW(I,K,J),K=1,14)
23  CONTINUE
    !END OF PROGRAM
100 CONTINUE
!print *, 'hi'

!FORMATS
101 FORMAT(1H ,/,20X,'WING LIFT DISTRIBUTION CALCULATION (WITH GROUND1 EFFECT)',/,20X,56('-'))
102 FORMAT(1H ,/,10X,'ALFA:',F10.2,8X,'LAMDA(1) :',F10.2,8X,'B :', F10.2,8X,'C:',F13.2,/,33X, &
'LAMDA(2) :',F10.2,8X,'S:', F10.2,8X,'AR:',F13.2,/,33X,'IB:',I10,8X,'JB :',I10,8X,'L.E. HEIGHT:', F6.2,/)
103 FORMAT(1H ,I3,' I ',F9.3,' II ',4(F9.3,' I '),' I ',4(F9.3,' I '))
104 FORMAT(1H ,'CL=',F10.4,2X,'L=',F10.4,4X,'CM=',F10.4,3X,'CD=',F10.4,3X,'L/L(INF)=',F10.4,4X,'GAMA/GAMA(INF)=',F10.4,/)
105 FORMAT(1H ,9X,'BB(',I3,')=',F10.4)
106 FORMAT(1H ,/,' T=',F10.2,3X,'SX=',F10.2,3X,'SZ=',F10.2,3X,'VT=',F10.2,3X,'TETA= ',F10.2,6X,'OMEGA=',F10.2)
107 FORMAT(1H ,//,' WAKE ELEMENTS,',//)
108 FORMAT(1H ,'QW(',I2,')=',22(F6.2))
109 FORMAT(1H ,' VORTIC(IT=',I3,')=',17(F6.3))
110 FORMAT(1H ,/,5X,'IDL',4X,'II',22X,'DCP',22X,'I I',25X,'GAMA',/,118('='),/,5X,'I',15X,'I= 1',11X,'2',11X,'3',&
11X,'4',5X,'I I',5X,'1',11X,'2',11X,'3',11X,'4',/,118('='))
111 FORMAT(1H ,9X,'ALF(',I2,')=',F10.4)
112 FORMAT(1H ,'QF(I=',I2,',J,X.Y.Z)= ',15(F6.1))
113 FORMAT(1H ,110('='))
STOP
END