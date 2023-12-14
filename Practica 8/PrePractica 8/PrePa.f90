PROGRAM P8

IMPLICIT NONE

DOUBLE PRECISION :: pi,x0,xf,phi0,phif,V,E,E1,E2,dphi0,dx,E3,E0
DOUBLE PRECISION,DIMENSION(2) :: funcin,dfuncout
DOUBLE PRECISION,DIMENSION(400) :: PHI,DPHI
INTEGER :: nequs,N,iteracio,i,j
COMMON/CONSTANTS/V,E
EXTERNAL EDO

pi = 4*datan(1.0d0)

x0 = 0.0d0
xf = 1.0d0
phi0 = 0.0d0
phif = 0.0d0
dphi0 = 0.15

V = -2.4

funcin = [phi0,dphi0]
N = 400
dx = xf/N
nequs = 2

open(15,file = 'resultatsprepra8.dat')


DO i =  1,4

E1 = (i*pi)**2/2.d0 + V +0.2d0
E2 = (i*pi)**2/2.d0 + V +0.3d0

CALL TIR(E1,E2,dx,funcin,nequs,EDO,N,E3,x0,xf,iteracio)
E0 = E3
CALL INTEGRAL_RK4(n,x0,xf,funcin,nequs,EDO,PHI,DPHI,E0)

WRITE(15,*) "#",E3

DO j =1,N
WRITE(15,*) x0+j*dx,PHI(j)
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

SUBROUTINE TIR(E1,E2,dx,funcin,nequs,edofuncio,N,E3,x_0,x_1,iteracio)

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(nequs) :: funcin,dfuncout1,dfuncout2,dfuncout3
DOUBLE PRECISION :: E1,E2,V,E,E3,dx,x_0,x_1,phif3,phif2,phif1
DOUBLE PRECISION,DIMENSION(N) :: PHI1,PHI2,PHI3,DPHI1,DPHI2,DPHI3
INTEGER :: nequs,N,i,iteracio
COMMON/CONSTANTS/V,E 
EXTERNAL :: edofuncio


DO i = 1,N

CALL INTEGRAL_RK4(n,x_0,x_1,funcin,nequs,edofuncio,PHI1,DPHI1,E1)
phif1 = PHI1(size(PHI1))
CALL INTEGRAL_RK4(n,x_0,x_1,funcin,nequs,edofuncio,PHI2,DPHI2,E2)
phif2 = PHI2(size(PHI2))


E3 = (E1*phif2-E2*phif1)/(phif2-phif1)

E=E3
CALL INTEGRAL_RK4(n,x_0,x_1,funcin,nequs,edofuncio,PHI3,DPHI3,E3)
phif3 = PHI3(size(PHI3))

IF (DABS(phif3).LT.(1.D-5)) THEN
E3 = E3
iteracio = i
EXIT
ELSE
E1 = E2
E2 = E3

END IF
END DO

END

SUBROUTINE INTEGRAL_RK4(n,x_0,x_1,funcin,nequs,edofuncio,PHI,DPHI,E0)

IMPLICIT NONE

INTEGER :: i,n,nequs,j
DOUBLE PRECISION :: dx,x_0,x_1,E0,V,E
DOUBLE PRECISION,DIMENSION(nequs) :: funcin,dfuncout
DOUBLE PRECISION,DIMENSION(N) :: PHI,DPHI
COMMON/CONSTANTS/V,E
EXTERNAL :: edofuncio

dx = (x_1-x_0)/n
E = E0
DO i = 1,n

CALL RungeKutta4order(dx,funcin,dfuncout,nequs,edofuncio)
PHI(i) = dfuncout(1)
DPHI(i) = dfuncout(2)
DO j = 1,nequs
funcin(j) = dfuncout(j)
END DO
END DO

END