PROGRAM P10

IMPLICIT NONE

DOUBLE PRECISION :: T_amb,T_Lt,T_0t,L,alpha,beta,tf,h,epsilon
DOUBLE PRECISION, allocatable :: theta1(:),theta2(:),theta3(:)
INTEGER :: nx,nt,i

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

allocate(theta1(nx),theta2(nx),theta3(nx))


h = L/dble(nx)
epsilon = 1d-6

OPEN(11,file="apartat1.dat")

theta1(:) = 0.d0
theta1(nx) = T_Lt - T_amb
beta = 0.00004

CALL RESOLUCIO(h,nx,theta1,epsilon,alpha,beta)
WRITE(11,*) "#"
do i =0,Nx
WRITE(11,*) i*h,theta1(i) + T_amb
END DO
WRITE(11,"(/)")

theta2(:) = 0.d0
theta2(nx) = T_Lt - T_amb
beta = 0.0003

CALL RESOLUCIO(h,nx,theta2,epsilon,alpha,beta)
WRITE(11,*) "#"
do i =0,Nx
WRITE(11,*) i*h,theta2(i) + T_amb
END DO
WRITE(11,"(/)")

theta3(:) = 0.d0
theta3(nx) = T_Lt - T_amb
beta = 0.0025

CALL RESOLUCIO(h,nx,theta3,epsilon,alpha,beta)
WRITE(11,*) "#"
do i =0,Nx
WRITE(11,*) i*h,theta3(i) + T_amb
END DO
WRITE(11,"(/)")

CLOSE(11)
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
WRITE(*,*) "Ha convergit en ",k,"iteracions"
EXIT
END IF

END DO

END
