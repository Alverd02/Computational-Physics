PROGRAM P7

IMPLICIT NONE

DOUBLE PRECISION :: pi,m,l,g,w_N,T_N,t_f,t_0,y_0,p_0,E_t_1,E_t_2,Ecine,Epoten
INTEGER :: n,i
DOUBLE PRECISION,DIMENSION(1500) :: theta_a_1,theta_a_2,p_a_1,p_a_2
DOUBLE PRECISION,DIMENSION(1500) :: theta_b_1,theta_b_2,p_b_1,p_b_2
DOUBLE PRECISION,DIMENSION(1500) :: theta_c_1,theta_c_2,p_c_1,p_c_2
COMMON/CONSTANTS/g,l,m
EXTERNAL edo

pi = 3.14159265359

m = 0.98
l = 1.07
g = 10.44
w_N = dsqrt(g/l)
T_N = (2.d0*pi)/w_N
t_0 = 0.d0
t_f = 7.d0*T_N

OPEN(11,file="P7-23-24-res.dat")
! 1)
n = 1500
y_0 = 0.025
p_0 = 0

CALL euler(n,p_0,y_0,t_0,t_f,theta_a_1,p_a_1,edo)
CALL euler_segon(n,p_0,y_0,t_0,t_f,theta_a_2,p_a_2,edo)

WRITE(11,*) "# 0)"
DO i = 1,n
WRITE(11,*) (t_0+(t_f-t_0)*i/n),theta_a_1(i),theta_a_2(i)
END DO

WRITE(11,"(/)")

n = 1500
y_0 = pi - 0.15
p_0 = 0

CALL euler(n,p_0,y_0,t_0,t_f,theta_b_1,p_b_1,edo)
CALL euler_segon(n,p_0,y_0,t_0,t_f,theta_b_2,p_b_2,edo)

WRITE(11,*) "# 1)"
DO i = 1,n
WRITE(11,*) (t_0+(t_f-t_0)*i/n),theta_b_1(i),p_b_1(i),theta_b_2(i),p_b_2(i)
END DO

WRITE(11,"(/)")

n = 1500
y_0 = pi - 0.025
p_0 = 0.12

CALL euler(n,p_0,y_0,t_0,t_f,theta_c_1,p_c_1,edo)
CALL euler_segon(n,p_0,y_0,t_0,t_f,theta_c_2,p_c_2,edo)

WRITE(11,*) "# 2)"

DO i = 1,n
E_t_1 = Ecine(p_c_1(i)) + Epoten(theta_c_1(i))
E_t_2 = Ecine(p_c_2(i)) + Epoten(theta_c_2(i))
WRITE(11,*) (t_0+(t_f-t_0)*i/n),E_t_1,Ecine(p_c_1(i)),E_t_2,Ecine(p_c_2(i))
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

SUBROUTINE euler_segon(n,p_0,y_0,x_0,x_f,y,p,funci)

IMPLICIT NONE

INTEGER :: n,i 
DOUBLE PRECISION :: p_0,y_0,x_0,x_f,funci,h,k_1,k_2
DOUBLE PRECISION,DIMENSION(n) :: y,p


h = (x_f-x_0)/n
y(1) = y_0
p(1) = p_0
DO i =2,n


y(i) = y(i-1) + (h/3.d0)*p(i-1) + (2*h/3.d0)*funci(y(i-1) + (3*h/4.d0)*p(i-1))
p(i) = p(i-1) + (h/3.d0)*funci(y(i-1)) + (2*h/3.d0)*funci(y(i-1) + (3*h/4.d0)*p(i-1))
END DO
END

DOUBLE PRECISION FUNCTION Ecine(x)

IMPLICIT NONE

DOUBLE PRECISION :: y,x,g,l,m
COMMON/CONSTANTS/g,l,m

y = (1/2.d0)*m*(x**2)*(l**2)
return
END

DOUBLE PRECISION FUNCTION Epoten(x)

IMPLICIT NONE

DOUBLE PRECISION :: y,x,g,l,m
COMMON/CONSTANTS/g,l,m

y = -m*g*l*dcos(x)
return
END