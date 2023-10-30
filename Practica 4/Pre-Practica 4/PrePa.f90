PROGRAM P4

IMPLICIT NONE

INTEGER :: i,h
DOUBLE PRECISION :: x_a_DOUBLE,x_b_DOUBLE
REAL :: x_a_single,x_b_single


! Apartat 0

! a)
OPEN(11,file="taules")
x_a_DOUBLE = 0.d0
WRITE(11,*) "Cas a): "
DO i=0,200000000


if ( mod (i, 100000).eq.0) THEN 
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
x = x+i*h
CALL func(x,valor)
resultat2 = valor + resultat2
END DO

resultat = resultat1 * h*resultat2
RETURN

END

SUBROUTINE  simpsontresvuit(x1,x2,k,func,resultat)

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,x,resultat2,resultat1,valor,valor_0,valor_N
INTEGER :: k,intervals,i

intervals = 3**k
h = (x2-x1)/intervals
x = x1

CALL func(x,valor_0)
CALL func(x,valor_N)
resultat1 = (h/3.)*(valor_0+valor_N)

resultat2 = 0.d0

DO i=1,(intervals-1)

x = x+i*h
CALL func(x,valor)
IF (MOD(i,2).EQ.0) THEN

resultat2 = resultat2 + 2.*valor

ELSE

resultat2 = resultat2 + 4.*valor

END IF

END DO

resultat = resultat1 + (h/3.)*resultat2

RETURN
END 