SUBROUTINE MATRX(A,N,G)
    !MATRX IS A MATRIX REDUCER OF THE GAUSSIAN TYPE
    !A(I,J) IS THE MATRIX, A(I,N) IS THE RHS VECTOR
    !AND G(I) IS THE SOLUTION VECTOR.
    REAL A(400,400),TEMP(400,400),G(400)
    !INITIALIZE THE G VECTOR TO ALL ZEROES
    DO I=1,N-1
    G(I)=0
    END DO
    !CONVERT COEFFICIENT MATRIX TO
    !UPPER TRIANGULAR FORM
    DO I=1,N-1
    5 IF(ABS(A(I,I)).LT.0.0000001) GOTO 9
    P=A(I,I)
    DO J=I,N
    A(I,J)=A(I,J)/P
    END DO
    DO K=I+1,N-1
    P2=A(K,I)
    DO L=I,N
    A(K,L)=A(K,L)-P2*A(I,L)
    END DO
    END DO
    END DO
    !BACK SUBSTITUTE TRIANGULARIZED MATRIX TO GET
    !VALUES OF SOLUTION VECTOR
    DO I=N-1,1,-1
    G(I)=A(I,N)
    DO J=1,N-1
    A(I,I)=0
    G(I)=G(I)-A(I,J)*G(J)
    END DO
    END DO
    RETURN
    !ORDER MATRIX SO THAT DIAGONAL COEFFICIENTS ARE
    !NOT =0 AND STOP IS MATRIX IS SINGULAR
    9 IF(I.NE.N-1) THEN
    DO J=1,N
    TEMP(I,J)=A(I,J)
    A(I,J)=A(I+1,J)
    A(I+1,J)=TEMP(I,J)
    END DO
    GOTO 5
    ELSE
    GOTO 10
    END IF
    10 WRITE(6,*) 'NO SOLUTION'
    STOP
    END


REAL EP(400,2),PT1(400,2),PT2(400,2),TH(400)
REAL CO(400,2),A(400,400),B(400,400),G(400)
REAL EPT(400,2),SIG(400),PHI(400),DL(400)
OPEN(8,FILE='CPLV.DAT',STATUS='REPLACE')
OPEN(11, FILE='A.txt', STATUS='REPLACE')
OPEN(12, FILE='G.txt', STATUS='REPLACE')
OPEN(13, FILE='CO.txt', STATUS='REPLACE')
OPEN(14, FILE='RHS.txt', STATUS='REPLACE')
OPEN(15, FILE='THETAS.txt', STATUS='REPLACE')
OPEN(16, FILE='B.txt', STATUS='REPLACE')
OPEN(17, FILE='sig.txt', STATUS='REPLACE')
OPEN(8,FILE='CPSD.DAT',STATUS='REPLACE')
OPEN(9,FILE='AFOIL2.DAT',STATUS='OLD')
WRITE(6,*) 'ENTER NUMBER OF PANELS'
READ(5,*) M
N=M+1
WRITE(6,*) 'ENTER ANGLE OF ATTACK IN DEGREES'
READ(5,*) ALPHA
AL=ALPHA/57.2958
!READ IN THE PANEL END POINTS
DO I=1,M+1
READ(9,*) EPT(I,1), EPT(I,2)
END DO
!CONVERT PANELING TO CLOCKWISE
DO I=1,M+1
EP(I,1)=EPT(N-I+1,1)
EP(I,2)=EPT(N-I+1,2)
END DO
!ESTABLISH COORDINATES OF PANEL END POINTS
DO I=1,M
PT1(I,1)=EP(I,1)
PT2(I,1)=EP(I+1,1)
PT1(I,2)=EP(I,2)
PT2(I,2)=EP(I+1,2)
END DO
!FIND PANEL ANGLES TH(J)
DO I=1,M
DZ=PT2(I,2)-PT1(I,2)
DX=PT2(I,1)-PT1(I,1)
TH(I)=ATAN2(DZ,DX)
END DO
!ESTABLISH SOURCE STRENGTHS (SIGMA=V DOT N)
DO I=1,M
SIG(I)=(COS(AL)*SIN(TH(I))-SIN(AL)*COS(TH(I)))
END DO
!ESTABLISH SURFACE POINTS (COLLOCATION POINTS)
DO I=1,M
CO(I,1)=(PT2(I,1)-PT1(I,1))/2+PT1(I,1)
CO(I,2)=(PT2(I,2)-PT1(I,2))/2+PT1(I,2)
END DO

! DEBUG:
!Write Thetas
DO I=1,M
    WRITE(15,*) TH(I)
END DO
! Write colocation points
do, i=1,m
    WRITE(13,*) CO(I,1),' ,',CO(I,2)
enddo
!Write source strengths
DO I=1,M
    WRITE(17,*) SIG(I)
END DO



!ESTABLISH INFLUENCE COEFFICIENTS
DO I=1,M
TEMP=0
DO J=1,M
!C CONVERT THE COLLOCATION POINT TO LOCAL PANEL COORDS.
XT=CO(I,1)-PT1(J,1)
ZT=CO(I,2)-PT1(J,2)
X2T=PT2(J,1)-PT1(J,1)
Z2T=PT2(J,2)-PT1(J,2)
X=XT*COS(TH(J))+ZT*SIN(TH(J))
Z=-XT*SIN(TH(J))+ZT*COS(TH(J))
X2=X2T*COS(TH(J))+Z2T*SIN(TH(J))
Z2=0
!SAVE PANEL LENGTHS
IF(I.EQ.1) THEN
DL(J)=X2
END IF
!COMPUTE R AND THETA VALUES FOR THE COLOC. POINT
R1=SQRT(X**2+Z**2)
R2=SQRT((X-X2)**2+Z**2)
TH1=ATAN2(Z,X)
TH2=ATAN2(Z,X-X2)

!COMPUTE THE DOUBLET INFLUENCE COEFFICIENTS
IF(I.EQ.J) THEN
A(I,J)=0.5
ELSE
A(I,J)=-0.15916*(TH2-TH1)
END IF


!COMPUTE THE SOURCE INFLUENCE COEFF'S AND ADD THEM UP TO GIVE THE RHS
IF(I.EQ.J) THEN
TEMP=TEMP+SIG(J)/3.14159265*(X*LOG(R1))
ELSE
TEMP=TEMP+SIG(J)/6.28319*(X*LOG(R1) * -(X-X2)*LOG(R2)+Z*(TH2-TH1))
END IF
END DO

!ADD WAKE INFLUENCE COEFF.
XW=CO(I,1)-PT2(M,1)
ZW=CO(I,2)-PT2(M,2)
DTHW=-ATAN(ZW/XW)
A(I,N)=-0.15916*(DTHW)

A(I,N+1)=TEMP
END DO
!ADD AN EXPLICIT KUTTA CONDITION
DO I=1,N+1
A(N,I)=0
END DO
A(N,1)=-1
A(N,M)=1
A(N,N)=-1

! Write A matrix
DO I=1,N
    DO J=1,N + 1
       write(11, '(F16.10)', advance='no') A(I,J)
    end do
    write(11, *) ''  ! this gives you the line break
 end do

 ! Write B matrix
DO I=1,M
   DO J=1,N
      write(16, '(F16.10)', advance='no') B(I,J)
   end do
   write(16, *) ''  ! this gives you the line break
end do



 ! WRITE RHS vector
 DO I=1,N
    WRITE(14, '(F16.10)') A(I,N+1)
 END DO

! SOLVE FOR THE SOLUTION VECTOR OF DOUBLET STRENGTHS

N=N+1
CALL MATRX(A,N,G)

!CONVERT DOUBLET STRENGTHS INTO TANGENTIAL
!VELOCITIES ALONG THE AIRFOIL SURFACE AND CP'S
!ON EACH PANEL.

200 CONTINUE
DO I=1,M
PHI(I)=CO(I,1)*COS(AL)+CO(I,2)*SIN(AL)+G(I)
END DO
DO I=1,M-1
R=(DL(I+1)+DL(I))/2
VEL=(PHI(I)-PHI(I+1))/R
CP=1-VEL**2
WRITE(8,*) PT2(I,1),', ',CP
END DO
WRITE(6,*) ' '
WRITE(6,*) 'LIFT COEFFICIENT=',G(M+1)
STOP

END