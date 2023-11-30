PROGRAM P7

IMPLICIT NONE

DOUBLE PRECISION :: m,l,g,w_N,t_0,t_f,T_N,pi
INTEGER :: n
COMMON/CONSTANTS/g,l
EXTERNAL edo

pi = 3.14159265359

m = 510.d0
l = 45.d0
g = 3.71
w_N = dsqrt(g/l)
T_N = (2.d0*pi)/w_N
t_0 = 0.d0
t_f = 6.d0*T_N
p_0 = 0

n = 1300

OPEN()

END PROGRAM P7

SUBROUTINE euler(n,p_0,y_0,x_0,x_f,y,funci)

IMPLICIT NONE

INTEGER :: n
DOUBLE PRECISION,DIMENSION(n) :: y,dfunci,p
DOUBLE PRECISION :: p_0,y_0,y_0,funci,x_0,x_f,h

CALL derivataula_edo(n,x_0,x_f,funci,dfunci)

h = (x_1-x_0)/n

y(1) = y_0
p(1) = p_0

DO i = 2,n

p(i) = p(i-1) + h*dfunci(i-1)
y(i) = y(i-1) + h*p(i-1)

END DO

END

SUBROUTINE derivataula_edo(ndates,x_0,x_f,funci,dfunci)
IMPLICIT NONE
INTEGER :: ndates,i
DOUBLE PRECISION :: h,funci,x_0,x_1,valorx
DOUBLE PRECISION, DIMENSION(1:ndates) :: dfunci

h = (x_1-x_0)/ndates

DO i=1,ndates
valorx = x_0 + i*h
dfunci(i) = funci(valorx)
END DO

END DO

DOUBLE PRECISION FUNCTION edo(y,x)

IMPLICIT NONE

DOUBLE PRECISION :: y,x,g,l
COMMON/CONSTANTS/g,l

y = (-g/l)*dsin(x)

END