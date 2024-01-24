! Interpolació

SUBROUTINE interpol(tin,xout)

COMMON /DADES/ times, positions
DOUBLE PRECISION, DIMENSION(1:501) :: times,positions
DOUBLE PRECISION :: xout,tin
INTEGER :: i

DO i=0,size(times)
if (times(i) <= tin .and. tin <= times(i + 1)) then
xout = positions(i) + (tin-times(i)) / (times(i+1) - times(i)) * (positions(i+1) - positions(i))
EXIT
END IF
END DO
RETURN 
END



! Trobar 0 de la funció fa falta el valor de la funció i la seva derivada

SUBROUTINE newtonraphson(fun,xini,preci,nitera,valorarrel)
IMPLICIT NONE
DOUBLE PRECISION :: xini,preci,valorarrel,x1,valordf,valorf,eps,t
INTEGER :: nitera,i


DO i=1,999999

CALL fun(xini,valorf,eps,valordf,t)   
x1 = xini - valorf/valordf
IF (ABS(valorf/valordf).LE.preci) THEN
valorarrel = x1
nitera = i

RETURN 
STOP
END IF
valorarrel = x1
nitera = i
xini = x1
END DO
END

! Metode de biseccio per trobar 0  de funcions, nomes va falta el valor de la funció

SUBROUTINE biseccio(fun,A,B,preci,nitera,valorarrel)
IMPLICIT NONE
DOUBLE PRECISION :: A,B,preci,valorarrel,valorf_A,valorf_B,d,valordf_A,valordf_B,valorf_d,valordf_d
INTEGER :: i,max_nitera,nitera

max_nitera=INT(LOG((B-A)/preci)/LOG(2.0d0))+1

IF (valorf_A*valorf_B.GT.0.0d0) THEN
WRITE(*,*) "No se puede aplicar el metodo, no hay cambio de singo"
STOP
END IF

DO i=1,max_nitera
d = (A+B)/2.0d0
CALL fun(A,valorf_A,valordf_A)
CALL fun(B,valorf_B,valordf_B)
CALL fun(d,valorf_d,valordf_d)

IF (valorf_d.EQ.0.0d0) THEN
valorarrel = d
nitera = i
RETURN 
STOP
END IF

IF ((valorf_A*valorf_d).LT.0.0D0) THEN
B = d
ELSE 
A = d
END IF


IF ((B-A).LT.preci) THEN
valorarrel = d
nitera = i
RETURN 
STOP
END IF

END DO

END

! Derivada primer ordre, numerica

SUBROUTINE derivataula(ndates,valorsx,funci,dfunci)
IMPLICIT NONE
INTEGER :: ndates,i
DOUBLE PRECISION :: h
DOUBLE PRECISION, DIMENSION(1:ndates) :: valorsx,funci,dfunci

h = valorsx(2)-valorsx(1)

dfunci(1) = (funci(2)-funci(1))/h

dfunci(ndates) = (funci(ndates)-funci(ndates-1))/h

DO i=2,ndates-1
dfunci(i) = (funci(i+1)-funci(i-1))/(2*h)

END DO

END 


! Meteode de trapezis

SUBROUTINE trapezoidalrule(x1,x2,k,func,resultat)

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,x,valor,resultat2,resultat1,valor_0,valor_N,func
INTEGER :: k,intervals,i

intervals = 3.d0**k
h = (x2-x1)/intervals

x=x1
valor_0 = func(x)
x = x2
valor_N =  func(x)

resultat1 = ((valor_0+valor_N)*h)/2.
resultat2 = 0.d0

DO i=1,(intervals-1)
x = x1
x = x+i*h
valor =  func(x)
resultat2 = valor + resultat2
END DO

resultat = resultat1 + h*resultat2
RETURN

END

! Calcula el valor d'una integral, d'una funció (retorna un valor), tambe es por arreglar perque la funció d'entrada sigui una llita de valors.


SUBROUTINE  simpsontresvuit(x1,x2,k,func,resultat)

! Metode de Simpsin 3/8 repetit, utilitzem 3^k intervals i começant des de 1 fins intervals-1 
!multipliquem per 2  o per 3 segons si i es dsivisible per 3 o no.

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,valor,x,valor0,valorN,func
INTEGER :: intervals,k,i

intervals = 3.d0**k
h = (x2-x1)/intervals

valor0 = func(x1)
valorN = func(x2)

resultat = valor0 + valorN

DO i=1,intervals-1
x = x1
IF (MOD(i,3).eq.0) THEN

valor =  func(x+i*h)
resultat=resultat+2*valor

ELSE

valor =  func(x+i*h)
resultat=resultat+3*valor

END IF

END DO

resultat = (3*h*resultat)/8.

RETURN
END 

!Genera uns histograma

subroutine histograma(ndat, xi, a, b, nbox, xhis, vhis, errhis, boxsize, ierr)
	implicit none
	integer :: ndat, nbox, ierr, i, ibox, icount
	double precision :: a, b, boxsize
	double precision, dimension(ndat) :: xi
	double precision, dimension(nbox) :: xhis, vhis, errhis

	if (a.ge.b) then 
		ierr=1
		return
	endif
	boxsize=(b-a)/nbox

	icount=0

	do i=1,nbox
		vhis(i)=0
		errhis(i)=0
	enddo

	do i=1,ndat
		if (xi(i).ge.a.and.xi(i).le.b) then 
			ibox=int((xi(i)-a)/boxsize)+1
		if (ibox.eq.nbox+1) ibox=nbox 

			vhis(ibox)=vhis(ibox)+1
			icount=icount+1
		endif
	enddo

	if (icount.eq.0) then 
		ierr=2
		return
	endif

	ierr=0
	print*,"accepted:",icount," out of:",ndat

	do i=1,nbox
		xhis(i)=a+boxsize/2.d0+(i-1)*boxsize
		errhis(i)=sqrt(vhis(i)/icount*(1.d0-vhis(i)/icount))/boxsize / sqrt(dble(icount))
		vhis(i)=vhis(i)/icount/boxsize
	enddo
end subroutine histograma

! Genera numeros aleatoris segons una distribució

SUBROUTINE accepta(ndades,numeros,xlow,xhigh,cotasup,funcio)

IMPLICIT NONE

DOUBLE PRECISION :: xlow,xhigh,cotasup,funcio,p
DOUBLE PRECISION, dimension(ndades) :: numeros
INTEGER :: ndades,n,i

n = 1

DO WHILE (n.lt.ndades)
numeros(n) = (xhigh-xlow)*RAND() + xlow
p = cotasup*RAND()
IF (funcio(numeros(n)).ge.p) then

n = n + 1

END IF
END DO

RETURN
END

! Genera numeros aleatoris gaussians

SUBROUTINE boxmuller(ndades,sigma,mu,yres)

IMPLICIT NONE
COMMON/CONSTANTS/pi
DOUBLE PRECISION,dimension(ndades) :: yres
INTEGER :: ndades,i
DOUBLE PRECISION :: mu,sigma,pi,x1,x2

DO i=1,ndades
x1 = RAND()
x2 = RAND()
yres(i) = mu + sigma*dsqrt(-2*dlog(x1))*dcos(2*pi*x2)
END DO

END

! Montecarlo cru vol dir que li pasem una funcio, els seus limits i retorna un valor i un error

SUBROUTINE montecarlocru(a,b,n,func,resultat,error)

IMPLICIT NONE

INTEGER :: n,i
DOUBLE PRECISION :: h,a,b,func,resultat,error,x2

resultat = 0.d0
x2 = 0.d0

DO i=1,n 

h = (b-a)*func((b-a)*rand()+a)
resultat = resultat + h
x2 = x2 + h**2

END DO

error = (1/dsqrt(dble(n)))*dsqrt((x2/dble(n)) - (resultat/dble(n))**2)
resultat = resultat/dble(n)

RETURN
END

! Integral de montecarlo: si la distribucio de "valors" es d'una desitibució dins de f(x) el meteode es queda com està, si f(x) es general i tenim valors d'una distribució externa
! h se l'haura de dividir per la distribució en el valor que toqui.( retorna un valor i el seru error)

subroutine montecarlo(valors,n,func,resultat,error)

IMPLICIT NONE

INTEGER :: n,i
DOUBLE PRECISION :: h,func,resultat,error,x2
DOUBLE PRECISION, DIMENSION(n) :: valors

resultat = 0.d0
x2 = 0.d0

DO i=1,n 

h = func(valors(i))
resultat = resultat + h
x2 = x2 + h**2

END DO

error = (1/dsqrt(dble(n)))*dsqrt((x2/dble(n)) - (resultat/dble(n))**2)
resultat = resultat/dble(n)

RETURN
END

! MonteCarlo multidimesnional

subroutine montecarlomulti(valors,n,func,resultat,error)

IMPLICIT NONE

INTEGER :: n,i,variables
DOUBLE PRECISION :: h,func,resultat,error,x2
DOUBLE PRECISION, DIMENSION(n) :: valors

resultat = 0.d0
x2 = 0.d0

DO i=5,n,5 

h = func(valors(i-1),valors(i-2),valors(i-3),valors(i-4),valors(i))
resultat = resultat + h
x2 = x2 + h**2

END DO

error = (1/dsqrt(dble(n)))*dsqrt((x2/dble(n)) - (resultat/dble(n))**2)
resultat = resultat/dble(n)

RETURN
END

DOUBLE PRECISION FUNCTION func6(x1,x2,x3,x4,x5)

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,x3,x4,x5,y,pi
COMMON/CONSTANTS/pi

y = exp(x2*dcos(x2+x3-x5))*(pi*x3**2*x4**2*x5**3+(dcos(x3+x4-2*x1))**2*x3*dsin(x5))
RETURN
END

! Metode de Euler per resoldre equs primer ordre

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

! Igual pero un metode millor

SUBROUTINE segonordre(T, x0, xf, y0, f0, y, dy, f)

IMPLICIT NONE

DOUBLE PRECISION :: x0, y0, f0, h, xf, f, dk1, k2, dk2, k1
DOUBLE PRECISION :: l, g
INTEGER :: T, i 
DOUBLE PRECISION, DIMENSION(T+1) :: y, dy


y(1) = y0 
dy(1) = f0

h = abs(xf-x0)/dble(T)

DO i=2, T+1

	k1 = dy(i-1)
    dk1 = f(y(i-1))
    k2 = dy(i-1) + 3.d0*h*dk1/4.d0
    dk2 = f(y(i-1)+ 3.d0*h*k1/4.d0)
	y(i) = y(i-1) + h*k1/3.d0 + 2.d0*h*k2/3.d0 
    dy(i) =  dy(i-1)  + h*dk1/3.d0  + 2.d0*h*dk2/3.d0
ENDDO

RETURN
END SUBROUTINE segonordre

! Runge-Kutta 4 resol un pas d'una equació diferencial d'un ordre de magnitud igual o mes que 2,

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

DOUBLE PRECISION :: 
DOUBLE PRECISION ,DIMENSION(nequs)::funcin,dyoutput
INTEGER :: nequs
 

dyoutput(1) = funcin(2)
dyoutput(2) =  ! equació diferecial


END

SUBROUTINE integralRK4(x0,xf, N, E0, phi, phi_f)
IMPLICIT NONE
DOUBLE PRECISION :: x0,xf,E0,dx,x,E,V0,phi_f,delta,const
DOUBLE PRECISION,DIMENSION(2) :: funcin,phiRK
DOUBLE PRECISION,DIMENSION(N) ::  phi
INTEGER :: N, nequs, i
COMMON/CONSTANTS/E,V0,delta,const
EXTERNAL EDO
funcin = [0.0d0,2*1.0D-6]
nequs = 2
E = E0
dx = (xf-x0)/dble(N) 

DO i =1,n
x = x0 + dx*i
phi(i) = funcin(1)
call RungeKutta4(x,dx,funcin,phiRK,nequs,EDO)
funcin = phiRK
END DO
phi_f = phiRK(1)
END

! Metode de tir

SUBROUTINE tir(E1, E2, N, error, E3)
IMPLICIT NONE
DOUBLE PRECISION :: E1, E2, E3, error, x0, xf, phi_f1, phi_f2, phi_f3,L
DOUBLE PRECISION,DIMENSION(N) :: phi
INTEGER :: N, i,iteracions

L = 14.0d0
x0 = -L/2.0d0
xf = L/2.0d0

iteracions = 10

DO i = 1, iteracions

CALL integralRK4(x0,xf, N, E1, phi, phi_f1)
CALL integralRK4(x0,xf, N, E2, phi, phi_f2)

E3 = (E1*phi_f2-E2*phi_f1)/(phi_f2-phi_f1)
        
CALL integralRK4(x0,xf, N, E3, phi, phi_f3)
WRITE(13,*) i,E3
IF (abs(phi_f3) .GT. error) then
E1 = E2
E2 = E3

ELSE
EXIT
END IF
END DO

END

! Resolucio de l'equació de Poisson en 2D, metode de Jacobi i Sobre-relax

SUBROUTINE RESOLUCIO(h,Nx,Ny,T,icontrol,epsilon,RHO)

IMPLICIT NONE

INTEGER :: i,j,Nx,Ny,k,icontrol
DOUBLE PRECISION :: T(Nx,Ny),RHO,epsilon,error,TNEW(Nx,Ny),TOLD(Nx,Ny),w,delta,h
TOLD = T
TNEW = TOLD

DO k = 1,100000


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

! Condicions de contorn d'una malla 2D condicions de Dirichlet, temperatures

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

! Solucio equs parabolica

SUBROUTINE Crank_Nicolson(theta0,alpha,beta,nx,dx,dt,nt,icontrol,T_amb)

IMPLICIT NONE

DOUBLE PRECISION :: theta0(0:nx),theta1(0:nx)
DOUBLE PRECISION :: r,beta,alpha,dx,dt,T_amb
DOUBLE PRECISION :: BB(1:nx-1,1:nx-1),AP(1:nx-1),A0(1:nx-1),Btheta(1:nx-1)
DOUBLE PRECISION :: AM(1:nx-1),thetax(1:nx-1)
INTEGER :: nx,i,nt,icontrol,k,j,n1,n2,n3,m,n,nmax,l

theta1 = theta0

BB = 0.d0
AP = 0.d0
A0 = 0.d0
AM = 0.d0

r = (alpha*dt)/dx**2

DO i = 1,nx-1


A0(i) = 2.*(1 + r) + beta*dt
BB(i,i) = 2.*(1 - r) - beta*dt

AP(i) = -r
AM(i) = -r

IF (i.NE.nx-1) THEN
BB(i,i+1) = r
BB(i+1,i) = r
END IF

END DO

AP(nx-1) = 0.d0
AM(1) = 0.d0

DO i=1,nt

DO j = 1,nx-1

Btheta(j) = 0.d0

DO k = 1,nx-1

Btheta(j) = Btheta(j) + BB(j,k)*theta0(k)

END DO
END DO

Btheta(1) = Btheta(1) + 2*r*theta1(0)
Btheta(nx-1) = Btheta(nx-1) + 2*r*theta1(nx)

DO n = 1, nx-1 
thetax(n) = theta0(n)
END DO
nmax = nx-1
CALL TRIDIAG(AM,A0,AP,Btheta,thetax,nmax)

DO m = 1, nx-1
theta0(m) = thetax(m)
END DO

IF (icontrol.EQ.1) THEN

n1 = int(0.06d0/dx)
n2 = int(0.42d0/dx)
n3 = int(1.26d0/dx)

WRITE(12,*) i*dt, theta0(n1) + T_amb, theta0(n2) + T_amb, theta0(n3) + T_amb

END IF

IF (icontrol.EQ.2) THEN

WRITE(13,*) i*dt, sum(theta0)/dble(nx) + T_amb

END IF

IF (icontrol.EQ.3)THEN
IF(mod(i, 10).EQ.0) THEN
DO l = 0,nx
WRITE(14,*) dt*i, dx*l, theta0(l) + T_amb
END DO
WRITE(14,"(/)")
END IF
END IF

END DO

END


SUBROUTINE TRIDIAG(A,B,C,R,PSI,IMAX)
IMPLICIT double precision (A-H,K,O-Z)
IMPLICIT INTEGER (I-J , L-N)
double precision  BET
double precision  GAM(4001)
double precision A(IMAX),B(IMAX),C(IMAX),R(IMAX),PSI(IMAX)

IF(B(1).EQ.0.) PAUSE
BET=B(1)
PSI(1)=R(1)/BET
DO 11 J=2,IMAX
GAM(J)=C(J-1)/BET
BET=B(J)-A(J)*GAM(J)
IF(BET.EQ.0) PAUSE
PSI(J)=(R(J)-A(J)*PSI(J-1))/BET
11      CONTINUE

DO 12 J=IMAX-1,1,-1
PSI(J)=PSI(J)-GAM(J+1)*PSI(J+1)
12      CONTINUE

RETURN
END
