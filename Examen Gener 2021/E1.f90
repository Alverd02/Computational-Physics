PROGRAM E1

IMPLICIT NONE
DOUBLE PRECISION :: M,L,g,a,b,w0,mu,resultat_monte,resultat_trap,error_monte,dtheta0,theta0,T_trap,T_monte,pi
INTEGER :: N,i
COMMON/CONSTANTS/M,L,g,w0,pi
COMMON/PARAMETRES/mu
EXTERNAL K

M = 1.25d0
L = 0.43d0
g = 9.8d0
w0 = dsqrt(g/L)
pi = 4.d0*datan(1.d0)

a = 0.d0
b = pi/2.d0

! 1
! a)

N = 350
dtheta0 = (pi-0.1d0)/45

OPEN(11,file="apartat1a.dat")
WRITE(11,*) "# Periode amb trapezis i montecarlo"
DO i = 1,45

mu = dsin(i*dtheta0/2.d0)**2

CALL trapezoidalrule(a,b,N,K,resultat_trap)
CALL montecarlocru(a,b,N,K,resultat_monte,error_monte)
T_trap = (4.d0/w0)*resultat_trap
T_monte = (4.d0/w0)*resultat_monte

WRITE(11,*) i*dtheta0,T_trap,T_monte,error_monte

END DO
WRITE(11,"(/)") 

! b)

WRITE(11,*) "# Convergencia"
theta0 = pi - 0.2d0
DO i=4,100,2

CALL trapezoidalrule(a,b,i,K,resultat_trap)
CALL montecarlocru(a,b,i,K,resultat_monte,error_monte)

T_trap = (4.d0/w0)*resultat_trap
T_monte = (4.d0/w0)*resultat_monte

WRITE(11,*) i,T_trap,T_monte,error_monte
END DO

CLOSE(11)

END PROGRAM E1

SUBROUTINE trapezoidalrule(x1,x2,k,func,resultat)

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,x,valor,resultat2,resultat1,valor_0,valor_N,func
INTEGER :: k,intervals,i

intervals = k
h = (x2-x1)/intervals

x=x1
valor_0 = func(x)
x = x2
valor_N =  func(x)

resultat1 = ((valor_0+valor_N)*h)/2.
resultat2 = 0.d0

DO i=1,(intervals-1)
x = x1
x = x+i*h
valor =  func(x)
resultat2 = valor + resultat2
END DO

resultat = resultat1 + h*resultat2
RETURN

END

SUBROUTINE montecarlocru(a,b,n,func,resultat,error)

IMPLICIT NONE

INTEGER :: n,i
DOUBLE PRECISION :: h,a,b,func,resultat,error,x2

resultat = 0.d0
x2 = 0.d0

DO i=1,n 

h = (b-a)*func((b-a)*rand()+a)
resultat = resultat + h
x2 = x2 + h**2

END DO

error = (1/dsqrt(dble(n)))*dsqrt((x2/dble(n)) - (resultat/dble(n))**2)
resultat = resultat/dble(n)

RETURN
END

DOUBLE PRECISION FUNCTION K(theta)

IMPLICIT NONE

DOUBLE PRECISION :: K,theta,mu
COMMON/PARAMETRES/mu
K = 1.d0/(dsqrt(1.d0-mu*dsin(theta)**2))

END