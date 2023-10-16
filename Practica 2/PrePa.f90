PROGRAM P1

IMPLICIT NONE

DOUBLE PRECISION :: radimano,R_k,L,t
INTEGER :: i,k

L = 25.d0

DO i=1,4

R_k = radimano(i,L)
write(*,*) R_k

END DO

END PROGRAM P1

DOUBLE PRECISION FUNCTION radimano(k,L)

IMPLICIT NONE

DOUBLE PRECISION :: L,R_k
INTEGER :: k

R_k = L - 0.15 - 0.3*(k-1)

RETURN

END 

SUBROUTINE posipisto(t,x)

IMPLICIT NONE

DOUBLE PRECISION :: t,r,omega,omega0,L,radimano
DOUBLE PRECISION :: DIMENSION(1:4) :: x
INTEGER :: k

L = 25.d0
omega0 = 4.8

DO k=1,4
r = radimano(k,L)
omega = omega0
END DO

END