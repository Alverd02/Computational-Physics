PROGRAM E2

IMPLICIT NONE
DOUBLE PRECISION :: M,L,g,a,b,w0,alpha,theta0,T,dtheta0_mes,dtheta0_menys,pi
DOUBLE PRECISION :: funcin(2)
DOUBLE PRECISION,ALLOCATABLE :: phi_mes(:),dphi_mes(:),phi_menys(:),dphi_menys(:),Ecin_menys(:)&
,Epot_menys(:),E_menys(:),Ecin_mes(:),Epot_mes(:),E_mes(:)
INTEGER :: N,i
COMMON/CONSTANTS/ M,L,g,w0,pi,alpha
EXTERNAL EDO

M = 1.25d0
L = 0.43d0
g = 9.8d0
w0 = dsqrt(g/L)
pi = 4.d0*datan(1.d0)
alpha = -0.1d0

T = 2.d0*pi*dsqrt(g/L)
theta0 = 0.d0
dtheta0_mes = 2*dsqrt(g/L) + 0.3d0
dtheta0_menys = 2*dsqrt(g/L) - 0.3d0

OPEN(11,file = "apartat2a.dat")

N = 5000

allocate(phi_mes(N),dphi_mes(N),phi_menys(N),dphi_menys(N))

funcin = [theta0,dtheta0_mes]
CALL integralRK4(funcin,0.d0,5.d0*T, N, phi_mes, dphi_mes)

funcin = [theta0,dtheta0_menys]
CALL integralRK4(funcin,0.d0,5.d0*T, N, phi_menys, dphi_menys)
WRITE(11,*)  "# Comparació"
DO i = 1,n

WRITE(11,*)  phi_mes(i),dphi_mes(i),phi_menys(i),dphi_menys(i)

END DO
WRITE(11,"(/)") 

allocate(Ecin_menys(N),Epot_menys(n),E_menys(n),Ecin_mes(N),Epot_mes(n),E_mes(n))

DO i=1,N


Ecin_menys(i) = M*L**2.d0*dphi_menys(i)**2/2.d0
Epot_menys(i) = M*L*g*(1-dcos(phi_menys(i)))
E_menys(i) = Ecin_menys(i)+Epot_menys(i)

Ecin_mes(i) = M*L**2.d0*dphi_mes(i)**2/2.d0
Epot_mes(i) = M*L*g*(1-dcos(phi_mes(i)))
E_mes(i) = Ecin_mes(i)+Epot_mes(i)

WRITE(11,*) (5*T/N)*i,Ecin_menys(i),Epot_menys(i),E_menys(i),Ecin_mes(i),Epot_mes(i),E_mes(i)

END DO
WRITE(11,"(/)") 

alpha = 0.d0

funcin = [theta0,dtheta0_mes]
CALL integralRK4(funcin,0.d0,5.d0*T, N, phi_mes, dphi_mes)
funcin = [theta0,dtheta0_menys]
CALL integralRK4(funcin,0.d0,5.d0*T, N, phi_menys, dphi_menys)

DO i=1,N


Ecin_menys(i) = M*L**2.d0*dphi_menys(i)**2/2.d0
Epot_menys(i) = M*L*g*(1-dcos(phi_menys(i)))
E_menys(i) = Ecin_menys(i)+Epot_menys(i)

Ecin_mes(i) = M*L**2.d0*dphi_mes(i)**2/2.d0
Epot_mes(i) = M*L*g*(1-dcos(phi_mes(i)))
E_mes(i) = Ecin_mes(i)+Epot_mes(i)

WRITE(11,*) (5*T/N)*i,Ecin_menys(i),Epot_menys(i),E_menys(i),Ecin_mes(i),Epot_mes(i),E_mes(i)

END DO

CLOSE(11)

END PROGRAM E2

SUBROUTINE RungeKutta4(x0,dx,funcin,dfuncout,nequs,edofuncio)
IMPLICIT NONE
DOUBLE PRECISION :: dx,x,x0
DOUBLE PRECISION,DIMENSION(nequs) :: funcin, dfuncout, y, k1, k2, k3, k4
INTEGER :: nequs



CALL edofuncio(nequs, x0, funcin, k1) !
y = funcin + 0.5d0 * k1 * dx 
x = x0 + 0.5d0*dx 
CALL edofuncio(nequs, x, y, k2)
y = funcin + 0.5d0 * k2 * dx
x = x0 +  0.5d0*dx
CALL edofuncio(nequs, x, y, k3)
y = funcin + k3 * dx 
x = x0 + dx 
CALL edofuncio(nequs, x, y, k4)

dfuncout = funcin + dx/6.d0*(k1+2.d0*k2+2.d0*k3+k4)

END

SUBROUTINE EDO(nequs, x,funcin , dyoutput)
IMPLICIT NONE

DOUBLE PRECISION :: M,L,g,w0,pi,alpha,x
DOUBLE PRECISION ,DIMENSION(nequs)::funcin,dyoutput
INTEGER :: nequs
COMMON/CONSTANTS/ M,L,g,w0,pi,alpha

dyoutput(1) = funcin(2)
dyoutput(2) =  -(g/L)*dsin(funcin(1)) + alpha*funcin(2)


END

SUBROUTINE integralRK4(funcin,x0,xf, N,phi,dphi)
IMPLICIT NONE
DOUBLE PRECISION :: x0,xf,dx,x
DOUBLE PRECISION,DIMENSION(2) :: funcin,phiRK
DOUBLE PRECISION,DIMENSION(N) ::  phi,dphi
INTEGER :: N, nequs, i
EXTERNAL EDO
nequs = 2
dx = (xf-x0)/dble(N) 

DO i =1,n
x = x0 + dx*i
phi(i) = funcin(1)
dphi(i) = funcin(2)
call RungeKutta4(x,dx,funcin,phiRK,nequs,EDO)
funcin = phiRK
END DO

END