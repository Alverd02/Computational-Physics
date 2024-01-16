PROGRAM P9

IMPLICIT NONE

DOUBLE PRECISION :: Lx,Ly,T0y,Txf,Tfy,Tx0,h,epsilon,Tint(3),x,y,RHO,rho1,rho2
INTEGER :: Nx,Ny,i,j,icontrol
DOUBLE PRECISION, ALLOCATABLE :: T(:,:)
EXTERNAL RHO,rho1,rho2

Lx = 18.5d0
Ly = 31.0d0

T0y = 2.2d0
Txf = 13.5d0
Tfy = 35.0d0
Tx0 = 4.16d0

h = 0.5d0

Nx = Lx/h+1
Ny = Ly/h+1
allocate(T(Nx,Ny))


epsilon = 1.D-6

Tint = [2.0d0,14.0d0,230.0d0]

OPEN(11,file="res.dat")

DO i = 1,3

icontrol = 0
CALL CONTORN(Nx,Ny,T0y,Tx0,Tfy,Txf,T,Tint(i))
WRITE(11,*) "# Tint Jacobi",Tint(i)
CALL RESOLUCIO(h,Nx,Ny,T,icontrol,epsilon)
WRITE(11,"(/)")


icontrol = 1
CALL CONTORN(Nx,Ny,T0y,Tx0,Tfy,Txf,T,Tint(i))
WRITE(11,*) "# Tint Sobrerelax",Tint(i)
CALL RESOLUCIO(h,Nx,Ny,T,icontrol,epsilon)
WRITE(11,"(/)")

END DO
CLOSE(11)

OPEN(12,file="mapa.dat")

icontrol = 1
CALL CONTORN(Nx,Ny,T0y,Tx0,Tfy,Txf,T,Tint(1))
CALL RESOLUCIO(h,Nx,Ny,T,icontrol,epsilon)
do i=1,Nx
x=dble(i)*h
do j=1,Ny
y=dble(j)*h
write(12,*) x,y,T(i,j)
enddo
write(12,*) ""

enddo

CLOSE(12)
END PROGRAM P9

double precision function rho1(x,y)
    implicit none
    double precision :: x, y, r


    r = sqrt((x-9.5d0)**2 + (y-12.0d0)**2)
    rho1 = 0.8*exp(-(r-2.d0)**2/0.35d0**2)

end function


double precision function rho2(x,y)
    implicit none
    double precision:: x, y
    
    rho2 = 0.d0

    !Fogó localitzat en un rectangle
    if ((x > 13) .and. (x < 15)) then 
        if ((y > 22) .and. (y < 26)) then
            rho2 = 1.1
        end if
    end if

end function


double precision function rho(x,y)
    implicit none
    double precision :: x,y, rho1, rho2

    rho = rho1(x,y) + rho2(x,y)

end function

SUBROUTINE RESOLUCIO(h,Nx,Ny,T,icontrol,epsilon)

IMPLICIT NONE

INTEGER :: i,j,Nx,Ny,k,icontrol
DOUBLE PRECISION :: T(Nx,Ny),RHO,epsilon,error,TNEW(Nx,Ny),TOLD(Nx,Ny),w,delta,h
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
w = 1.52
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

WRITE(11,*) k,TNEW(16,27)

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
do i=1,nx
 T(i,1)=Tx0
 T(i,ny)=Txf
enddo
do j=1,ny
 T(1,j)=T0y
 T(nx,j)=Tfy
enddo
! inicialitzem les temperatures interiors a 1
do i=2,nx-1
 do j=2,ny-1
 T(i,j)=Tint
 enddo
enddo

END