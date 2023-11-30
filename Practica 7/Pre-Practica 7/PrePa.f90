PROGRAM P7

IMPLICIT NONE

DOUBLE PRECISION :: m,l,g,w_N,t_0,t_f,T_N,pi,p_0,y_0,Ecine,Epoten,E_P_1,E_P_2,E_T_1,E_T_2,E_C_1,E_C_2
DOUBLE PRECISION,DIMENSION(1300) :: p_a_1,theta_a_1,p_a_2,theta_a_2
DOUBLE PRECISION,DIMENSION(1800) :: p_b,theta_b
DOUBLE PRECISION,DIMENSION(2500) :: p_c_1,theta_c_1,p_c_2,theta_c_2
INTEGER :: n,i
COMMON/CONSTANTS/g,l,m
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

CALL euler(n,p_0,y_0,t_0,t_f,theta_a_1,p_a_1,edo)
CALL Adams_Bashforth(n,p_0,y_0,t_0,t_f,theta_a_2,p_a_2,edo)
WRITE(11,*) "# 0)"
DO i = 1,n
WRITE(11,*) (t_0+(t_f-t_0)*i/n),p_a_1(i),p_a_2(i)
END DO

WRITE(11,"(/)")

n = 1800
y_0 = pi-0.025
p_0 = 0

CALL euler(n,p_0,y_0,t_0,t_f,theta_b,p_b,edo)
WRITE(11,*) "# 1)"
DO i = 1,n
WRITE(11,*) (t_0+(t_f-t_0)*i/n),p_b(i),theta_b(i)
END DO
WRITE(11,"(/)") 

n = 2500
y_0 = 1.d0
p_0 = 0

CALL euler(n,p_0,y_0,t_0,t_f,theta_c_1,p_c_1,edo)

y_0 = pi-0.042

CALL euler(n,p_0,y_0,t_0,t_f,theta_c_2,p_c_2,edo)
WRITE(11,*) "# 2)"
DO i=1,n 

E_P_1 = Epoten(p_c_1(i))
E_C_1 = Ecine(theta_c_1(i))

E_T_1 = E_P_1 + E_C_1

E_P_2 = Epoten(p_c_2(i))
E_C_2 = Ecine(theta_c_2(i))

E_T_2 = E_P_2 + E_C_2

WRITE(11,*) (t_0+(t_f-t_0)*i/n),E_P_1,E_T_1,E_P_2,E_T_2

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

DOUBLE PRECISION FUNCTION Ecine(x)

IMPLICIT NONE

DOUBLE PRECISION :: y,x,g,l,m
COMMON/CONSTANTS/g,l,m

y = (1/2.d0)*m*x**2*l**2
return
END

DOUBLE PRECISION FUNCTION Epoten(x)

IMPLICIT NONE

DOUBLE PRECISION :: y,x,g,l,m
COMMON/CONSTANTS/g,l,m

y = -m*g*l*dcos(x)
return
END

SUBROUTINE Adams_Bashforth(n,p_0,y_0,x_0,x_f,y,p,funci)

IMPLICIT NONE

INTEGER :: n,i 
DOUBLE PRECISION :: p_0,y_0,x_0,x_f,funci,h
DOUBLE PRECISION,DIMENSION(n) :: y,p

CALL EULER(n,p_0,y_0,x_0,x_f,y,p,funci)

h = (x_f-x_0)/n

DO i =1,n-2 

y(i+2) = y(i+1) - (1/2.d0)*h*p(i) + (3/2.d0)*h*p(i+1)
p(i+2) = p(i+1) - (1/2.d0)*h*funci(y(i)) + (3/2.d0)*h*funci(y(i+1))
END DO

END