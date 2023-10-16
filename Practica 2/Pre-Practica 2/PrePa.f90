PROGRAM P1

IMPLICIT NONE

DOUBLE PRECISION :: radimano,L,r,t
DOUBLE PRECISION, DIMENSION(1:4) :: x
INTEGER :: k,i

! Apartat 1

L = 25.d0

DO k=1,4
r = radimano(k,L)
write(*,*) "Radi",k,":",r
END DO

! Apartat 2
t = 0.1
CALL posipisto(t,x)
write(*,*) "Vector:",x

! Apartat 3
OPEN(11,file="P2-23-24-res1.dat")

DO i =0,500
t=0.01*i
CALL posipisto(t,x)
write(11,"(5F10.2)") t,x
END DO

CLOSE(11)

END PROGRAM P1

! Apartat 1

DOUBLE PRECISION FUNCTION radimano(k,L)

IMPLICIT NONE

INTEGER :: k
DOUBLE PRECISION :: L,r

r = L - 0.15 - 0.3*(k-1)

RETURN 
END 

! Apartat 2

SUBROUTINE posipisto(t,x)

DOUBLE PRECISION :: t,L,omega0,omega,radimano,r
DOUBLE PRECISION, DIMENSION(1:4) :: x
INTEGER :: k

L = 25.d0
omega0 = 4.8

DO k=1,4
omega = omega0*(k/3.25+1)
r = radimano(K,L)
x(k) = r*cos(omega*t) + sqrt(L**2-r**2*(sin(omega*t))**2)
END DO 

RETURN 
END