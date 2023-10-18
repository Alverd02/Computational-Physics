PROGRAM P1

IMPLICIT NONE

DOUBLE PRECISION :: radimano,L,r,t,x1,x2,x3,x4,xout
DOUBLE PRECISION, DIMENSION(1:4) :: x
DOUBLE PRECISION, DIMENSION(1:501) :: times,positions
INTEGER :: k,i,iostat
COMMON /DADES/ times, positions

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

! Apartat 6

OPEN(11,file="P2-23-24-res1.dat")
DO i=1,9999
READ(11,*,IOSTAT=iostat) t,x1,x2,x3,x4
IF(iostat<0) THEN
EXIT
ELSE
times(i) = t
positions(i) = x2
END IF
END DO
CLOSE(11)

! Apartat 7

OPEN(12,file="P2-23-24-res2.dat")
xout = 0
DO i =0,1200
t=0.0025*i
CALL interpol(t,xout)
write(12,"(5F10.2)") t,xout
END DO

WRITE(*,*) size(times)

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

! Apartat 6


SUBROUTINE interpol(tin,xout)

COMMON /DADES/ times, positions
DOUBLE PRECISION, DIMENSION(1:501) :: times,positions
DOUBLE PRECISION :: xout,tin
INTEGER :: i

DO i=0,size(times)
if (times(i) <= tin .and. tin <= times(i + 1)) then
xout = positions(i) + (tin-times(i)) / (times(i+1) - times(i)) * (positions(i+1) - positions(i))
EXIT
END IF
END DO
RETURN 
END