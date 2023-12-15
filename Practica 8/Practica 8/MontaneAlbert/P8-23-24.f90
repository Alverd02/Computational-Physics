PROGRAM P8

IMPLICIT NONE

DOUBLE PRECISION :: E,V0,delta,const,L,x0,xf,E1,E2,E3,E4,phi_f,error,E0,E5,E6
DOUBLE PRECISION,DIMENSION(400) :: PHI1,PHI2,PHI3,PHI4,PHI1_N,PHI2_N,PHI3_N,PHI4_N
INTEGER :: n,i,iteracio
COMMON/CONSTANTS/E,V0,delta,const
EXTERNAL EDO

V0 = -50.0d0
delta = 0.4d0
const = 3.80995
L = 14.0d0

x0 = -L/2.0d0
xf = L/2.0d0
n = 400


OPEN(12,file="integral.dat")
! 1)
E1 = -31.0D0
CALL integralRK4(x0,xf, N, E1, PHI1, phi_f)
WRITE(12,*) "# E = ",E1
DO i = 1,n

WRITE(12,*) x0 + ((xf-x0)*i)/n,PHI1(i)
END DO

WRITE(12,"(/)")

E2 = -30.0D0
CALL integralRK4(x0,xf, N, E2, PHI2, phi_f)
WRITE(12,*) "# E = ",E2
DO i = 1,n

WRITE(12,*) (x0 + ((xf-x0)*i)/n),PHI2(i)
END DO
WRITE(12,"(/)")

E3 = -14.0D0
CALL integralRK4(x0,xf, N, E3, PHI3, phi_f)
WRITE(12,*) "# E = ",E3
DO i = 1,n

WRITE(12,*) (x0 + ((xf-x0)*i)/n),PHI3(i)
END DO
WRITE(12,"(/)")

E4 = -13.0D0
CALL integralRK4(x0,xf, N, E4, PHI4, phi_f)
WRITE(12,*) "# E = ",E4
DO i = 1,n

WRITE(12,*) (x0 + ((xf-x0)*i)/n),PHI4(i)
END DO
WRITE(12,"(/)")


CLOSE(12)

OPEN(13,file="convergencia.dat")

! 2)

error = 1.d-6
WRITE(13,*) "#"
CALL tir(E1, E2, N, error, E0)
WRITE(13,"(/)")
WRITE(13,*) "#"
CALL tir(E3, E4, N, error, E0)
WRITE(13,"(/)")
WRITE(13,*) "#"
E5 = -4.0D0
E6 = -3.5D0
CALL tir(E5, E6, N, error, E0)

CLOSE(13)

! 2)

OPEN(14,file="integral_norm.dat")

E1 = -31.0D0
CALL integralRK4(x0,xf, N, E1, PHI1, phi_f)
CALL normalitzar(X0,XF,N,PHI1,PHI1_N)
WRITE(14,*) "# E = ",E1
DO i = 1,n

WRITE(14,*) x0 + ((xf-x0)*i)/n,PHI1_N(i)
END DO

WRITE(14,"(/)")

E2 = -30.0D0
CALL integralRK4(x0,xf, N, E2, PHI2, phi_f)
CALL normalitzar(x0,xf, N,PHI2,PHI2_N)
WRITE(14,*) "# E = ",E2
DO i = 1,n

WRITE(14,*) (x0 + ((xf-x0)*i)/n),PHI2_N(i)
END DO
WRITE(14,"(/)")

E3 = -14.0D0
CALL integralRK4(x0,xf, N, E3, PHI3, phi_f)
CALL normalitzar(x0,xf, N,PHI3,PHI3_N)
WRITE(14,*) "# E = ",E3
DO i = 1,n

WRITE(14,*) (x0 + ((xf-x0)*i)/n),PHI3_N(i)
END DO
WRITE(14,"(/)")

E4 = -13.0D0
CALL integralRK4(x0,xf, N, E4, PHI4, phi_f)
CALL normalitzar(x0,xf, N,PHI4,PHI4_N)
WRITE(14,*) "# E = ",E4
DO i = 1,n

WRITE(14,*) (x0 + ((xf-x0)*i)/n),PHI4_N(i)
END DO
WRITE(14,"(/)")


CLOSE(14)

END

SUBROUTINE RungeKutta4(x0,dx,funcin,dfuncout,nequs,edofuncio)
IMPLICIT NONE
DOUBLE PRECISION :: dx,x,x0
DOUBLE PRECISION,DIMENSION(nequs) :: funcin, dfuncout, y, k1, k2, k3, k4
INTEGER :: nequs
CALL edofuncio(nequs, x0, funcin, k1) !
y = funcin + 0.5d0 * k1 * dx 
x = x0 + 0.5d0*dx 
CALL edofuncio(nequs, x, y, k2)
y = funcin + 0.5d0 * k2 * dx
x = x0 +  0.5d0*dx
CALL edofuncio(nequs, x, y, k3)
y = funcin + k3 * dx 
x = x0 + dx 
CALL edofuncio(nequs, x, y, k4)
dfuncout = funcin + dx/6.d0*(k1+2.d0*k2+2.d0*k3+k4)

END

SUBROUTINE EDO(nequs, x,funcin , dyoutput)
IMPLICIT NONE
DOUBLE PRECISION :: x, E, V,V0,delta,const
DOUBLE PRECISION ,DIMENSION(nequs)::funcin,dyoutput
INTEGER :: nequs
COMMON/CONSTANTS/E,V0,delta,const


V = (V0*dsinh(2.0d0))/(dcosh(2.0d0)+dcosh(x/delta))    

dyoutput(1) = funcin(2)
dyoutput(2) =  ((V-E) * funcin(1))/const

END

SUBROUTINE integralRK4(x0,xf, N, E0, phi, phi_f)
IMPLICIT NONE
DOUBLE PRECISION :: x0,xf,E0,dx,x,E,V0,phi_f,delta,const
DOUBLE PRECISION,DIMENSION(2) :: funcin,phiRK
DOUBLE PRECISION,DIMENSION(N) ::  phi
INTEGER :: N, nequs, i
COMMON/CONSTANTS/E,V0,delta,const
EXTERNAL EDO
funcin = [0.0d0,2*1.0D-6]
nequs = 2
E = E0
dx = (xf-x0)/dble(N) 

DO i =1,n
x = x0 + dx*i
phi(i) = funcin(1)
call RungeKutta4(x,dx,funcin,phiRK,nequs,EDO)
funcin = phiRK
END DO
phi_f = phiRK(1)
END

SUBROUTINE tir(E1, E2, N, error, E3)
IMPLICIT NONE
DOUBLE PRECISION :: E1, E2, E3, error, x0, xf, phi_f1, phi_f2, phi_f3,L
DOUBLE PRECISION,DIMENSION(N) :: phi
INTEGER :: N, i,iteracions

L = 14.0d0
x0 = -L/2.0d0
xf = L/2.0d0

iteracions = 10

DO i = 1, iteracions

CALL integralRK4(x0,xf, N, E1, phi, phi_f1)
CALL integralRK4(x0,xf, N, E2, phi, phi_f2)

E3 = (E1*phi_f2-E2*phi_f1)/(phi_f2-phi_f1)
        
CALL integralRK4(x0,xf, N, E3, phi, phi_f3)
WRITE(13,*) i,E3
IF (abs(phi_f3) .GT. error) then
E1 = E2
E2 = E3

ELSE
EXIT
END IF
END DO

END

sUBROUTINE  simpson38(x1,x2,k,func,resultat)

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,valor,x,valor0,valorN
double precision,dimension(k) :: func
INTEGER :: intervals,k,i

intervals = k
h = (x2-x1)/intervals

VALOR0 = func(1)
VALORn = func(k)

resultat = valor0 + valorN

DO i=1,intervals-1
x = x1
IF (MOD(i,3).eq.0) THEN

valor = func(i+1)
resultat=resultat+2*valor

ELSE

valor =  func(i+1)
resultat=resultat+3*valor

END IF

END DO

resultat = (3*h*resultat)/8.

RETURN
END 

SUBROUTINE normalitzar(a,b, N,PHI,PHI_N)
IMPLICIT NONE
DOUBLE PRECISION,DIMENSION(N) :: phi2, PHI_n,phi
DOUBLE PRECISION :: x, dx, phi_F,resultat, a, b
INTEGER :: N, i

dx = (b-a)/n
DO i = 1,N
phi2(i) = phi(I)**2
END DO
    

CALL simpson38(a, b, N, phi2, resultat)

DO i = 1,N
PHI_n(i) = PHI(i)/sqrt(resultat)
END DO
END