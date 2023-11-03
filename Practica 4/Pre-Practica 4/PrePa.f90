PROGRAM P4

IMPLICIT NONE

INTEGER :: i,h,k
DOUBLE PRECISION :: x_a_DOUBLE,x_b_DOUBLE,pi,e,resultat_t_a,resultat_s_a,resultat_t_b,resultat_s_b,L,h2,error_t,error_s,t,resultat_t
DOUBLE PRECISION :: resultat_s,resultat_s4,resultat_t4,resultat_4s_b,resultat_4t_b
REAL :: x_a_single,x_b_single,intervals
COMMON/CONSTANTS/pi,e
EXTERNAL corba
EXTERNAL densitat_lineal
EXTERNAL f3

pi = 3.1415926535898
e = 2.718281828459
! Apartat 0

! a)
OPEN(11,file="taules")
x_a_DOUBLE = 0.d0
WRITE(11,*) "Cas a): "
DO i=0,200000000


IF ( MOD(i,100000).eq.0) THEN 
WRITE(11,"(e14.8)") x_a_DOUBLE
END IF 

x_a_DOUBLE = x_a_DOUBLE + 0.02

END DO

! b)

WRITE(11,*) "Cas b): "
h = 1000

DO i=0,2000

x_b_DOUBLE = 2*h*i
WRITE(11,"(e14.8)") x_b_DOUBLE

END DO
CLOSE(11)

! Haurien de ser iguals, pero veiem que hi ha discrepancies ja que la estrategia b) està arrodonint. 
!Per tant, fent servir la estrategia b) estem perdent xifres significatives, serà més adient la estrategia a).

! Apartat 2

OPEN(12,file="P4-23-24-res1.dat")

k = 13

CALL trapezoidalrule(-pi,pi,k,corba,resultat_t_a)
CALL simpsontresvuit(-pi,pi,k,corba,resultat_s_a)

WRITE(12,*) "Area sota la corba"
WRITE(12,*) "Metode de trapezis", resultat_t_a
WRITE(12,*) "Metode de simpson 3/8", resultat_s_a
WRITE(12,"(/)") 
L = 17.76

CALL simpsontresvuit(-L,L,k,densitat_lineal,resultat_s_b)
CALL trapezoidalrule(-L,L,k,densitat_lineal,resultat_t_b)

WRITE(12,*) "Massa de la barra"
WRITE(12,*) "Metode de trapezis", resultat_t_b
WRITE(12,*) "Metode de simpson 3/8", resultat_s_b

CLOSE(12)

! Apartat 3

OPEN(13,file="P4-23-24-res2.dat")

k = 13

CALL trapezoidalrule(-pi,pi,k,corba,resultat_t)
CALL simpsontresvuit(-pi,pi,k,corba,resultat_s)

DO i=2,13

intervals = 3**i
h2 = (2*pi)/intervals

CALL trapezoidalrule(-pi,pi,i,corba,resultat_t_a)
CALL simpsontresvuit(-pi,pi,i,corba,resultat_s_a)

error_t = (2*pi)*h2**2
error_s = (2*pi)*h2**4

WRITE(13,*) h2,error_t,error_s,abs(resultat_t-resultat_t_a),abs(resultat_s-resultat_s_a)
END DO

CLOSE(13)

OPEN(14,file="P4-23-24-res3.dat")

k = 13

CALL trapezoidalrule(-pi,pi,k,densitat_lineal,resultat_t)
CALL simpsontresvuit(-pi,pi,k,densitat_lineal,resultat_s)

DO i=2,13

intervals = 3**i
L = 17.76
h2 = (2*L)/intervals


CALL trapezoidalrule(-pi,pi,i,densitat_lineal,resultat_t_b)
CALL simpsontresvuit(-pi,pi,i,densitat_lineal,resultat_s_b)

error_t = (2*L)*h2**2
error_s = (2*L)*h2**4

WRITE(14,*) h2,error_t,error_s,abs(resultat_t-resultat_t_b),abs(resultat_s-resultat_s_b)
END DO

CLOSE(14)

! Apartat 4

OPEN(15,file="P4-23-24-res4.dat")

k = 13

CALL trapezoidalrule(-pi/2.,pi/2.,k,densitat_lineal,resultat_t)
CALL simpsontresvuit(-pi/2.,pi/2,k,densitat_lineal,resultat_s)
CALL trapezoidalrule(-pi/2.,pi/2.,k,f3,resultat_t4)
CALL simpsontresvuit(-pi/2.,pi/2,k,f3,resultat_s4)

DO i=2,10

intervals = 3**i
h2 = pi/intervals

CALL trapezoidalrule(-pi/2.,pi/2.,i,densitat_lineal,resultat_t_b)
CALL simpsontresvuit(-pi/2.,pi/2,i,densitat_lineal,resultat_s_b)
CALL trapezoidalrule(-pi/2.,pi/2.,i,f3,resultat_4t_b)
CALL simpsontresvuit(-pi/2.,pi/2,i,f3,resultat_4s_b)

WRITE(15,*) h2,error_t,error_s,abs(resultat_t-resultat_t_b),abs(resultat_s-resultat_s_b),abs(resultat_t4-resultat_4t_b), &
            abs(resultat_s4-resultat_4s_b)
END DO

CLOSE(15)

END PROGRAM P4


! Apartat 1

SUBROUTINE trapezoidalrule(x1,x2,k,func,resultat)

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,x,valor,resultat2,resultat1,valor_0,valor_N
INTEGER :: k,intervals,i

intervals = 3**k
h = (x2-x1)/intervals

x=x1
CALL func(x,valor_0)
x = x2
CALL func(X,valor_N)

resultat1 = ((valor_0+valor_N)*h)/2.
resultat2 = 0.d0

DO i=1,(intervals-1)
x = x1
x = x+i*h
CALL func(x,valor)
resultat2 = valor + resultat2
END DO

resultat = resultat1 + h*resultat2
RETURN

END

SUBROUTINE  simpsontresvuit(x1,x2,k,func,resultat)

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,valor,x,valor0,valorN
INTEGER :: intervals,k,i

intervals = 3**k
h = (x2-x1)/intervals

CALL func(x1,valor0)
CALL func(x2,valorN)

resultat = valor0 + valorN

DO i=1,intervals-1
x = x1
IF (MOD(i,3).eq.0) THEN

CALL func(x+i*h,valor)
resultat=resultat+2*valor

ELSE

CALL func(x+i*h,valor)
resultat=resultat+3*valor

END IF

END DO

resultat = (3*h*resultat)/8.

RETURN
END 

! Apartat 2

SUBROUTINE corba(x,valor)

IMPLICIT NONE

COMMON/CONSTANTS/pi,e
DOUBLE PRECISION :: valor,x,A_0,e,pi

A_0 = 0.33
valor = A_0*(cos(x-2.)*e**(-x**2-sin(x)))**2*sqrt(pi-x)

RETURN
END

SUBROUTINE densitat_lineal(x,valor)

IMPLICIT NONE


DOUBLE PRECISION :: valor,x,ro,L

ro = 0.72
L = 17.76
valor = ro*sqrt(1-(x/L)**2)*(1-(x/L))*((x/L)**2 + (x/L) + 1)

RETURN
END

SUBROUTINE f3(t,valor)

IMPLICIT NONE


DOUBLE PRECISION :: valor,t,ro,L

ro = 0.72
L = 17.76
valor = ro*sqrt(1-sin(t)**2)*(1-sin(t))*(sin(t)**2 + sin(t) + 1)*L*cos(t)

RETURN
END