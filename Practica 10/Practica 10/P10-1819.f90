PROGRAM P10

IMPLICIT NONE

DOUBLE PRECISION :: T_amb,T_Lt,T_0t,L,alpha,beta,tf,dx,epsilon,dt
DOUBLE PRECISION, allocatable :: theta1(:),theta2(:),theta3(:)
INTEGER :: nx,nt,i,j

T_amb = 22.d0
T_Lt = 280.d0
T_0t = T_amb

L = 1.50
alpha = 2.2d-5
beta = 0.00017

nx = 50
tf = 4500.d0
nt = 6000

! 1) Estacionari 

allocate(theta1(0:nx),theta2(0:nx),theta3(0:nx))

dx = L/dble(nx)
epsilon = 1d-6

OPEN(11,file="apartat1.dat")

theta1(:) = 0.d0
theta1(nx) = T_Lt - T_amb
beta = 0.00004

CALL RESOLUCIO(dx,nx,theta1,epsilon,alpha,beta)
WRITE(11,*) "#"
do i =0,Nx
WRITE(11,*) i*dx,theta1(i) + T_amb
END DO
WRITE(11,"(/)")

theta2(:) = 0.d0
theta2(nx) = T_Lt - T_amb
beta = 0.0003

CALL RESOLUCIO(dx,nx,theta2,epsilon,alpha,beta)
WRITE(11,*) "#"
do i =0,Nx
WRITE(11,*) i*dx,theta2(i) + T_amb
END DO
WRITE(11,"(/)")

theta3(:) = 0.d0
theta3(nx) = T_Lt - T_amb
beta = 0.0025

CALL RESOLUCIO(dx,nx,theta3,epsilon,alpha,beta)
WRITE(11,*) "#"
do i =0,Nx
WRITE(11,*) i*dx,theta3(i) + T_amb
END DO
WRITE(11,"(/)")

CLOSE(11)

! 2)
! a)

dt = tf/nt


END PROGRAM P10

SUBROUTINE RESOLUCIO(h,nx,T,epsilon,alpha,beta)

IMPLICIT NONE

INTEGER :: i,Nx,k
DOUBLE PRECISION :: T(nx),epsilon,error,TOLD,alpha,beta,h
DO k = 1,999999
error = 0.0d0

DO i = 1,nx-1
Told = T(i)
T(i) = (T(i+1) + T(i-1))/(2 + (h**2)*beta/alpha)

IF (ABS(Told-T(i)).GT.error) THEN
error = ABS(Told-T(i))
END IF

END DO
IF (error.LT.epsilon) THEN
EXIT
END IF

END DO

END

SUBROUTINE Crank_Nicolson(theta0,alpha,beta,nx,dx,dt,nt,icontrol,T_amb)

IMPLICIT NONE

DOUBLE PRECISION :: theta0(0:nx),theta1(0:nx)
DOUBLE PRECISION :: r,beta,alpha,dx,dt,T_amb
DOUBLE PRECISION :: BB(1:nx-1,1:nx-1),AP(1:nx-1),A0(1:nx-1),Btheta(1:nx-1)
DOUBLE PRECISION :: AM(nx-1)
INTEGER :: nx,i,nt,icontrol,k,j,ni,n2,n3

theta1 = theta0


r = (alpha*dt)/dx**2

DO i = 1,nx-1


A0(i) = 2.*(1 + r) + beta*dt
BB(i,i) = 2.*(1 - r) + beta*dt

AP(i) = -r
AM(i) = -r

IF (i.NE.nx-1) THEN
BB(i,i+1) = r
BB(i+1,i) = r
END IF

END DO

AP(1) = 0.d0
AM(nx-1) = 0.d0

DO i=1,nt

DO j = 1,nx-1

Btheta(j) = 0.d0

DO k = 1,nx-1

Btheta(j) = Btheta(j) + BB(j,k)*theta1(j)

END DO
END DO

Btheta(1) = Btheta(1) + 2*r*theta0(0)
Btheta(nx-1) = Btheta(nx-1) + 2*r*theta0(nx)

CALL TRIDIAG(AM,A0,AP,Btheta,theta1,nx-1)

IF (icontrol.EQ.1) THEN

if (iwrite == 1) then 

n1 = int(0.06d0/dx)
n2 = int(0.42d0/dx)
n3 = int(1.26d0/dx)

write(12,*) i*dt, theta1(n1) + T_amb, theta1(n2) + T_amb, theta1(n3) + T_amb

end if

END IF

END DO

END


SUBROUTINE TRIDIAG(A,B,C,R,PSI,IMAX)
IMPLICIT double precision (A-H,K,O-Z)
IMPLICIT INTEGER (I-J , L-N)
double precision  BET
double precision  GAM(4001)
double precision A(IMAX),B(IMAX),C(IMAX),R(IMAX),PSI(IMAX)

IF(B(1).EQ.0.) STOP
BET=B(1)
PSI(1)=R(1)/BET
DO 11 J=2,IMAX
GAM(J)=C(J-1)/BET
BET=B(J)-A(J)*GAM(J)
IF(BET.EQ.0) STOP
PSI(J)=(R(J)-A(J)*PSI(J-1))/BET
11      CONTINUE

DO 12 J=IMAX-1,1,-1
PSI(J)=PSI(J)-GAM(J+1)*PSI(J+1)
12      CONTINUE

RETURN
END