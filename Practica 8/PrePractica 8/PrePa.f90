PROGRAM P8

IMPLICIT NONE

DOUBLE PRECISION :: pi,x0,xf,phi0,phif,v
EXTERNAL EDO

pi = 4*datan(1.0d0)

x0 = 0.0d0
xf = 1.0d0
phi0 = 0.0d0
phif = 0.0d0

v = -2.4

END PROGRAM P8

SUBROUTINE RungeKutta4order(dx,funcin,dfuncout,nequs,edofuncio)

IMPLICIT NONE

INTEGER :: nequs,i
DOUBLE PRECISION,DIMENSION(nequs) :: funcin,dfuncout
DOUBLE PRECISION,DIMENSION(nequs) :: K1, K2, K3, K4

CALL edofuncio(funcin,K1)
dfuncout = funcin + dx*K1/2.D0

CALL edofuncio(dfuncout,K2)
dfuncout = funcin + dx*K2/2.D0

CALL edofuncio(dfuncout,K3)
dfuncout = funcin + dx*K3

CALL edofuncio(dfuncout,K4)
dfuncout = funcin + (dx/6.D0)*(K1+2.D0*K2+2.D0*K3+K4)

END

SUBROUTINE EDO(nequs,x,yinput,dyoutput)

IMPLICIT NONE

INTEGER :: nequs
DOUBLE PRECISION,DIMENSION(nequs) :: yinput,dyoutput
DOUBLE PRECISION :: x

DO i = 1,nequs

dyoutput(1) = yinput(2)
dyoutput(2) = -2.d0*(E*yinput(1)-V*yinput(1))

END