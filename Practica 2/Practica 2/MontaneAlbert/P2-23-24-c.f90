PROGRAM P2

IMPLICIT NONE

INTEGER :: i,iostat
DOUBLE PRECISION, DIMENSION(1:5) :: x
DOUBLE PRECISION :: pi,r,radimano,L,phi,phi_i,t,xout,xout_0,x1,x2,x3,x4,xinterpol,xinterpol0,x5
DOUBLE PRECISION, DIMENSION(1:82):: TI,XI
COMMON/CONSTANTS/pi
COMMON/POSIS/XI,TI

pi = 3.14159265

! Apartat 1

L = 18.5

DO i=1,5

r = radimano(i,L)
phi_i = phi(i)

WRITE(*,*) "Radi manovella",i,"=",r
WRITE(*,*) "Fase inicial",i,"=",phi_i

END DO

t = 0.1
CALL posipisto(t,x)
write(*,*) "Vector:",x

! Apartat 3

OPEN(11,file="P2-23-24-res1-c.dat")

DO i =0,81
t=0.0987654321*i
CALL posipisto(t,x)
write(11,"(6F10.2)") t,x
END DO

CLOSE(11)

! Apartat 6

OPEN(11,file="P2-23-24-res1-c.dat")
DO i=1,9999
READ(11,*,IOSTAT=iostat) t,x1,x2,x3,x4,x5
IF(iostat<0) THEN
EXIT
ELSE
TI(i) = t
XI(i) = x4
END IF
END DO
CLOSE(11)

! Apartat 7

OPEN(12,file="P2-23-24-res2-c.dat")

DO i =0,2000
t=0.003*i

xout =  xinterpol(t)

xout_0 = xinterpol0(t)

write(12,"(5F10.2)") t,xout_0,xout
END DO
CLOSE(12)

END PROGRAM P2

! Apartat 1

DOUBLE PRECISION FUNCTION radimano(i,L)

IMPLICIT NONE

INTEGER :: i
DOUBLE PRECISION :: L,r

r = L/i - 0.5

RETURN 
END

DOUBLE PRECISION FUNCTION phi(i)

IMPLICIT NONE

INTEGER :: i
DOUBLE PRECISION :: phi_i,pi 
COMMON/CONSTANTS/pi

phi_i = (i/5.)**2*pi

RETURN 
END

! Apartat 2

SUBROUTINE posipisto(t,x)

DOUBLE PRECISION :: t,L,omega,radimano,r,phi_i,phi
DOUBLE PRECISION, DIMENSION(1:5) :: x
INTEGER :: i

L = 18.5
omega = 5

DO i=1,5
r = radimano(i,L)
phi_i = phi(i)
x(i) = r*cos(omega*t + phi_i) + sqrt(L**2-r**2*(sin(omega*t + phi_i))**2)
END DO 

RETURN 
END

! Apartat 6 

DOUBLE PRECISION FUNCTION xinterpol(t)

COMMON /POSIS/ XI,TI
DOUBLE PRECISION, DIMENSION(1:82) :: XI,TI
DOUBLE PRECISION :: t,xout
INTEGER :: i

DO i=0,size(TI)
if (TI(i) <= t .and. t <= TI(i + 1)) then
xout = XI(i) + (t-TI(i)) / (TI(i+1) - TI(i)) * (XI(i+1) - XI(i))
EXIT
END IF
END DO
RETURN 
END

DOUBLE PRECISION FUNCTION xinterpol0(t)

COMMON /POSIS/ XI,TI
DOUBLE PRECISION, DIMENSION(1:82) :: XI,TI
DOUBLE PRECISION :: t,xout_0
INTEGER :: i

DO i=0,size(TI)
if (TI(i) <= t .and. t <= TI(i + 1)) then
xout_0 = XI(i)

EXIT

END IF
END DO
RETURN 
END