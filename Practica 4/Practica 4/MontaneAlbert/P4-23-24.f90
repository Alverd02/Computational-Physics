PROGRAM P4

IMPLICIT NONE

DOUBLE PRECISION :: a,b,valor_exacte,pi,x1,x2,h,A_T,A_S,error_t,error_s,A_i,A_nexti,A_m
INTEGER :: i,intervals
COMMON/CONSTANTS/a,b,pi
EXTERNAL YKohoutek

b = 429.074
a = 508.633
pi = 3.1415926535898

valor_exacte = a*b*(3*dsqrt(3.d0)+2.d0*pi)/24.d0

OPEN(11,file="P4-23-24-res.dat")

WRITE(11,"(A11)") "# Apartat a"

DO i=2,14

x2 = -(7/2.)*a
x1 = -4*a
intervals = 3.d0**i
h = (x2-x1)/intervals

CALL trapezoidalrule(x1,x2,i,YKohoutek,A_T)
CALL simpsontresvuit(x1,x2,i,YKohoutek,A_S)

WRITE(11,"(3e21.14)") h,A_T,A_S

END DO

WRITE(11,"(/)")

WRITE(11,"(A11)") "# Apartat b"

DO i=2,13

intervals = 3.d0**i
h = (x2-x1)/intervals

CALL trapezoidalrule(x1,x2,i,YKohoutek,A_T)
CALL simpsontresvuit(x1,x2,i,YKohoutek,A_S)

error_t = h**2/35.d0
error_s = h**4/10000000.d0

WRITE(11,"(5e21.14)") h,error_t,error_s,abs(valor_exacte-A_T),abs(valor_exacte-A_S)
END DO

WRITE(11,"(/)")
WRITE(11,"(A11)") "# Apartat c"

DO i=2,13

x2 = -(7/2.d0)*a
x1 = -4*a
intervals = 3.d0**i
h = (x2-x1)/intervals

CALL  trapezoidalrule(x1,x2,i,YKohoutek,A_i)
CALL  trapezoidalrule(x1,x2,i+1,YKohoutek,A_nexti)
A_m = (9*A_nexti-A_i)/8.d0

WRITE(11,"(3e21.14)") h,A_m,abs(valor_exacte-A_m)

END DO

! Podem observar com el metode dels trapezis modificat, és equivalent al Simpson 3/8

CLOSE(11)

END PROGRAM P4

SUBROUTINE trapezoidalrule(x1,x2,k,func,resultat)

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,x,valor,resultat2,resultat1,valor_0,valor_N,func
INTEGER :: k,intervals,i

intervals = 3.d0**k
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

SUBROUTINE  simpsontresvuit(x1,x2,k,func,resultat)

! Metode de Simpsin 3/8 repetit, utilitzem 3^k intervals i começant des de 1 fins intervals-1 
!multipliquem per 2  o per 3 segons si i es dsivisible per 3 o no.

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,valor,x,valor0,valorN,func
INTEGER :: intervals,k,i

intervals = 3.d0**k
h = (x2-x1)/intervals

valor0 = func(x1)
valorN = func(x2)

resultat = valor0 + valorN

DO i=1,intervals-1
x = x1
IF (MOD(i,3).eq.0) THEN

valor =  func(x+i*h)
resultat=resultat+2*valor

ELSE

valor =  func(x+i*h)
resultat=resultat+3*valor

END IF

END DO

resultat = (3*h*resultat)/8.

RETURN
END 

DOUBLE PRECISION FUNCTION YKohoutek(x)

IMPLICIT NONE

COMMON/CONSTANTS/a,b
DOUBLE PRECISION :: x,y,b,a 



y = b*sqrt(1-((x+4*a)**2/a**2))

RETURN
END