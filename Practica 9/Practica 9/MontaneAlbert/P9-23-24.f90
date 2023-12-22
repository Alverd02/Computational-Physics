PROGRAM P9

IMPLICIT NONE

DOUBLE PRECISION :: Lx,Ly,T0y,Txf,Tfy,Tx0,h,epsilon,Tint(3),x,y,RHO,RHO01,RHO02,RHO03
INTEGER :: Nx,Ny,i,j,icontrol
DOUBLE PRECISION, ALLOCATABLE :: T(:,:)
COMMON/CONSTANTS/RHO01,RHO02,RHO03
EXTERNAL RHO

RHO01 = 10.d0
RHO02 = 7.0d0
RHO03 = 5.5d0

Lx = 45.5d0
Ly = 33.5d0

T0y = 0.5d0
Txf = 25.3d0
Tfy = 11.2d0
Tx0 = 17.0d0

h = 0.25d0

Nx = Lx/h+1
Ny = Ly/h+1
allocate(T(Nx,Ny))


epsilon = 0.01

Tint = [15.0d0,220.0d0,1280.0d0]

OPEN(11,file="res.dat")



DO i = 1,3

icontrol = 0
CALL CONTORN(Nx,Ny,T0y,Tx0,Tfy,Txf,T,Tint(i))
WRITE(11,*) "# Tint Jacobi",Tint(i)
CALL RESOLUCIO(h,Nx,Ny,T,icontrol,epsilon,RHO)
WRITE(11,"(/)")


icontrol = 1
CALL CONTORN(Nx,Ny,T0y,Tx0,Tfy,Txf,T,Tint(i))
WRITE(11,*) "# Tint Sobrerelax",Tint(i)
CALL RESOLUCIO(h,Nx,Ny,T,icontrol,epsilon,RHO)
WRITE(11,"(/)")

END DO

icontrol = 1
CALL CONTORN(Nx,Ny,T0y,Tx0,Tfy,Txf,T,Tint(1))
CALL RESOLUCIO(h,Nx,Ny,T,icontrol,epsilon,RHO)
WRITE(11,"(/)")
DO i=1,Nx
X=dble(i)*h
DO j=1,Ny
y=dble(j)*h
WRITE(11,*) x,y,T(i,j)
END DO
WRITE(11,*) ""
END DO

RHO01 = 0
RHO02 = 0
RHO03 = 0

icontrol = 1
WRITE(11,*) ""
CALL CONTORN(Nx,Ny,T0y,Tx0,Tfy,Txf,T,Tint(1))
CALL RESOLUCIO(h,Nx,Ny,T,icontrol,epsilon,RHO)
WRITE(11,"(/)")
DO i=1,Nx
X=dble(i)*h
DO j=1,Ny
y=dble(j)*h
WRITE(11,*) x,y,T(i,j)
END DO
WRITE(11,*) ""
END DO

CLOSE(11)


END PROGRAM P9


double precision function rho(x,y)
    implicit none
    double precision :: rho1, rho2, rho3, x, y, r1, r3,RHO01,RHO02,RHO03
    COMMON/CONSTANTS/RHO01,RHO02,RHO03
    !Primer fogonet
    r1 = sqrt((x-22.5d0)**2 + (y-8d0)**2)
    rho1 = RHO01 * exp(-(r1-4.d0)**2/0.7**2)

    !Segon fogó
    rho2 = 0.d0
    if ((x>29) .and. (x<35) .and. (y>18) .and. (y<22)) then
        rho2 = RHO02
    end if

    !Tercer fogó
    r3 = sqrt((x-10.5d0)**2 + (y-22d0)**2)
    rho3 = RHO03 * exp(-(r3-5.d0)**2/1.2**2)

    !sumem
    rho = rho1 + rho2 + rho3

end function

SUBROUTINE RESOLUCIO(h,Nx,Ny,T,icontrol,epsilon,RHO)

IMPLICIT NONE

INTEGER :: i,j,Nx,Ny,h,k,icontrol
DOUBLE PRECISION :: T(Nx,Ny),RHO,epsilon,error,TNEW(Nx,Ny),TOLD(Nx,Ny),w,delta
TOLD = T
TNEW = TOLD

DO k = 1,999999

error = 0.0d0

DO i = 2,Nx-1
DO j = 2,Ny-1

IF (icontrol.EQ.0) THEN
TNEW(i,j) = (TOLD(i+1,j) + TOLD(i-1,j) + TOLD(i,j+1) + TOLD(i,j-1) + h**2*RHO(i*h,j*h))/4

IF (ABS(TNEW(i,j)-TOLD(i,j))>error) THEN
error = ABS(TNEW(i,j)-TOLD(i,j))
END IF

ELSE
w = 1.45
delta = 0.25*(TNEW(i+1,j) + TNEW(i-1,j) + TNEW(i,j+1) + TNEW(i,j-1) - 4*TNEW(i,j) + h**2*RHO(i*h,j*h))
TNEW(i,j) = TNEW(i,j) + w*delta

IF (ABS(delta)>error) THEN
error = ABS(delta)
END IF

END IF
END DO
END DO

IF (icontrol.EQ.0) THEN

TOLD = TNEW

END IF

WRITE(11,*) k,TNEW(102,54)

IF (error<epsilon) THEN
T = TNEW
EXIT
END IF

END DO

END

SUBROUTINE CONTORN(Nx,Ny,T0y,Tx0,Tfy,Txf,T,Tint)

IMPLICIT NONE
INTEGER :: Nx,Ny,i,j
DOUBLE PRECISION :: T(Nx,Ny),T0y,Tx0,Tfy,Txf,Tint
DO i=1,nx
T(i,1)=Tx0
T(i,ny)=Txf
END DO
DO j=1,ny
T(1,j)=T0y
T(nx,j)=Tfy
END DO

DO i=2,nx-1
DO j=2,ny-1
T(i,j)=Tint
END DO
END DO

END