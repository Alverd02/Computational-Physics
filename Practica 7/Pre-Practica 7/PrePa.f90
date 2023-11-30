PROGRAM P7

IMPLICIT NONE

DOUBLE PRECISION :: m,l,g,w_N,t_0,t_f,T_N,pi,p_0,y_0
DOUBLE PRECISION,DIMENSION(1300) :: p,theta
INTEGER :: n,i
COMMON/CONSTANTS/g,l
EXTERNAL edo

pi = 3.14159265359

m = 510.d0
l = 0.45
g = 3.71
w_N = dsqrt(g/l)
T_N = (2.d0*pi)/w_N
t_0 = 0.d0
t_f = 6.d0*T_N




OPEN(11,file="P7-23-24-res.dat")

n = 1300
y_0 = 0.02
p_0 = 0

CALL euler(n,p_0,y_0,t_0,t_f,theta,p,edo)
WRITE(11,*) "# 0)"
DO i = 1,n
WRITE(11,*) (t_0+(t_f-t_0)*i/n),p(i),theta(i)
END DO

WRITE(11,"(/)")

n = 1800
y_0 = pi-0.025

CALL euler(n,p_0,y_0,t_0,t_f,theta,p,edo)
WRITE(11,*) "# 1)"
DO i = 1,n
WRITE(11,*) (t_0+(t_f-t_0)*i/n),p(i),theta(i)
END DO
WRITE(11,"(/)") 

CLOSE(11)

END PROGRAM P7

SUBROUTINE euler(n,p_0,y_0,x_0,x_f,y,p,funci)

IMPLICIT NONE

INTEGER :: n,i
DOUBLE PRECISION,DIMENSION(n) :: y,p
DOUBLE PRECISION :: p_0,y_0,funci,x_0,x_f,h

h = (x_f-x_0)/n

y(1) = y_0
p(1) = p_0

DO i = 2,n

p(i) = p(i-1) + h*funci(y(i-1))
y(i) = y(i-1) + h*p(i-1)

END DO

END

DOUBLE PRECISION FUNCTION edo(x)

IMPLICIT NONE

DOUBLE PRECISION :: y,x,g,l
COMMON/CONSTANTS/g,l

y = (-g/l)*dsin(x)
return
END