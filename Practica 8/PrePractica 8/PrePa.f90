PROGRAM P8

IMPLICIT NONE

DOUBLE PRECISION :: pi,x0,xf,phi0,phif,V,E,E1,E2,dphi0,dx,E3
DOUBLE PRECISION,DIMENSION(2) :: funcin
DOUBLE PRECISION,DIMENSION(400) :: PHI,PHI_N
INTEGER :: nequs,N,iteracio,i,j
COMMON/CONSTANTS/V,E
EXTERNAL EDO

pi = 4*datan(1.0d0)

x0 = 0.0d0
xf = 1.0d0
phi0 = 0.0d0
phif = 0.0d0
dphi0 = 0.15d0

V = -2.4

funcin = [phi0,dphi0]

nequs = 2

open(15,file = 'resultatsprepra8.dat')

N = 400
dx = xf/N

DO i =  1,4

E1 = (i*pi)**2/2.d0 + V + 0.1
E2 = (i*pi)**2/2.d0 + V + 0.2

CALL TIR(E1,E2,dx,funcin,nequs,EDO,N,E3,x0,xf,iteracio,PHI)
CALL NORM(x0,xf,PHI,PHI_N,N)

WRITE(15,*) "# Iteracio = ",iteracio,"E = ",E3
DO j = 1,N

WRITE(15,*) x0+dx*j,PHI(j),PHI_N(j)

END DO
WRITE(15,"(/)")
END DO

N = 20
dx = xf/N

DO i =  1,4

E1 = (i*pi)**2/2.d0 + V + 0.1
E2 = (i*pi)**2/2.d0 + V + 0.2

CALL TIR(E1,E2,dx,funcin,nequs,EDO,N,E3,x0,xf,iteracio,PHI)
CALL NORM(x0,xf,PHI,PHI_N,N)

WRITE(15,*) "# Iteracio = ",iteracio,"E = ",E3
DO j = 1,N

WRITE(15,*) x0+dx*j,PHI(j),PHI_N(j)

END DO
WRITE(15,"(/)")
END DO

CLOSE(15)

END PROGRAM P8

SUBROUTINE RungeKutta4order(dx,funcin,dfuncout,nequs,edofuncio)

IMPLICIT NONE

INTEGER :: nequs
DOUBLE PRECISION,DIMENSION(nequs) :: funcin,dfuncout
DOUBLE PRECISION,DIMENSION(nequs) :: K1, K2, K3, K4
DOUBLE PRECISION :: dx,x
EXTERNAL :: edofuncio

CALL edofuncio(nequs,x,funcin,K1)
dfuncout = funcin + dx*K1/2.D0

CALL edofuncio(nequs,x,dfuncout,K2)
dfuncout = funcin + dx*K2/2.D0

CALL edofuncio(nequs,x,dfuncout,K3)
dfuncout = funcin + dx*K3

CALL edofuncio(nequs,x,dfuncout,K4)
dfuncout = funcin + (dx/6.D0)*(K1+2.D0*K2+2.D0*K3+K4)

END

SUBROUTINE EDO(nequs,x,yinput,dyoutput)

IMPLICIT NONE

INTEGER :: nequs
DOUBLE PRECISION,DIMENSION(nequs) :: yinput,dyoutput
DOUBLE PRECISION :: x,V,E
COMMON/CONSTANTS/V,E 

dyoutput(1) = yinput(2)
dyoutput(2) = -2.d0*(E*yinput(1)-V*yinput(1))

END

SUBROUTINE TIR(E1,E2,dx,funcin,nequs,edofuncio,N,E3,x_0,x_1,iteracio,PHI)

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(nequs) :: funcin,funcin1,dfuncout1,dfuncout2,dfuncout3,funcin2,funcin3
DOUBLE PRECISION,DIMENSION(N) :: PHI
DOUBLE PRECISION :: E1,E2,V,E,E3,dx,x_0,x_1
INTEGER :: nequs,N,i,iteracio
COMMON/CONSTANTS/V,E 
EXTERNAL :: edofuncio


funcin1 = funcin
funcin2 = funcin
funcin3 = funcin

iteracio = 1

DO WHILE (iteracio.LT.N)

DO i = 1,N
E = E1
CALL RungeKutta4order(dx,funcin1,dfuncout1,nequs,edofuncio)
funcin1 = dfuncout1
E  = E2
CALL RungeKutta4order(dx,funcin2,dfuncout2,nequs,edofuncio)
funcin2 = dfuncout2
END DO


E3 = (E1*dfuncout2(1)-E2*dfuncout1(1))/(dfuncout2(1)-dfuncout1(1))

DO i = 1,N

E=E3
PHI(i) = funcin3(1)
CALL RungeKutta4order(dx,funcin3,dfuncout3,nequs,edofuncio)
funcin3 = dfuncout3
END DO

IF ((DABS(dfuncout3(1))).LT.(1.d-5)) THEN

E1 = E2
E2 = E3
iteracio = iteracio + 1
ELSE

EXIT

END IF

END DO
END

SUBROUTINE NORM(x1,x2,PHI,PHI_N,N)

IMPLICIT NONE

INTEGER :: N,i
DOUBLE PRECISION,DIMENSION(N) :: PHI_N,PHI2,PHI
DOUBLE PRECISION :: V,E,resultat,x1,x2
COMMON/CONSTANTS/V,E 
EXTERNAL :: edofuncio

DO i = 1,N

PHI2(i) = PHI(i)**2

END DO

CALL simpsontresvuit(x1,x2,n,PHI2,resultat)

DO i = 1,N

PHI_N(i) = PHI(i)/DSQRT(resultat)

END DO

END

SUBROUTINE  simpsontresvuit(x1,x2,k,func,resultat)

! Metode de Simpsin 3/8 repetit, utilitzem 3^k intervals i começant des de 1 fins intervals-1 
!multipliquem per 2  o per 3 segons si i es dsivisible per 3 o no.

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,valor,x,valor0,valorN
DOUBLE PRECISION,DIMENSION(k) :: func
INTEGER :: intervals,k,i

intervals = k
h = (x2-x1)/intervals

valor0 = func(1)
valorN = func(size(func))

resultat = valor0 + valorN

DO i=1,intervals-1
x = x1
IF (MOD(i,3).eq.0) THEN

valor =  func(i)
resultat=resultat+2*valor

ELSE

valor =  func(i)
resultat=resultat+3*valor

END IF

END DO

resultat = (3*h*resultat)/8.

RETURN
END 