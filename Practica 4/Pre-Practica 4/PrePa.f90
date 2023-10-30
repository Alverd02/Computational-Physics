PROGRAM P4

IMPLICIT NONE

INTEGER :: i,h
DOUBLE PRECISION :: x_a_DOUBLE,x_b_DOUBLE

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

INTEGER :: i,h
REAL :: x_a_REAL,x_b_REAL

! Apartat 0

! a)
OPEN(12,file="taules_2")
x_a_REAL = 0.
WRITE(12,*) "Cas a): "
DO i=0,200000000


if ( mod (i, 100000).eq.0) THEN 
WRITE(12,"(e14.8)") x_a_REAL
END IF 

x_a_REAL = x_a_REAL + 0.02

END DO

! b)

WRITE(12,*) "Cas b): "
h = 1000

DO i=0,2000

x_b_REAL = 2*h*i
WRITE(12,"(e14.8)") x_b_REAL

END DO
CLOSE(12)

END PROGRAM P4