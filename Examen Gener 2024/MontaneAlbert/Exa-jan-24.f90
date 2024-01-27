PROGRAM Exam

IMPLICIT NONE

DOUBLE PRECISION :: t0,tf,lambda,omega,phi0,z0,x1,x2,resultat_tancat,xlow,xhigh,cotasup,error,resultat
DOUBLE PRECISION,ALLOCATABLE :: funcina(:),za(:),phia(:),funcinb(:),zb(:),phib(:),numeros(:)
INTEGER :: n,nequs,i,k,ISEED,ndades
COMMON/CONSTANTS/lambda,omega
EXTERNAL EDO
EXTERNAL f
EXTERNAL g
EXTERNAL h

ISEED=20470586
CALL SRAND(ISEED)

t0 = 0.d0
tf = 10.d0
lambda = 2.5d0
omega = 1.d0

phi0 = 0.d0



n = 400
nequs = 2

allocate(funcina(nequs),za(n),phia(n),funcinb(nequs),zb(n),phib(n))


OPEN(11,file="Exa-jan-24-res1.dat")
! 1)
! a)
! Cas A)

z0 = 0.3d0

funcina = [z0,phi0]

CALL integralRK4(t0,tf,n,nequs,funcina, za, phia)
WRITE(11,*) "# Cas A"
DO i=1,n

WRITE(11,*) za(i),phia(i)

END DO
WRITE(11,"(/)") 

! Cas B)

z0 = 0.9d0

funcinb = [z0,phi0]

CALL integralRK4(t0,tf,n,nequs,funcinb, zb, phib)

WRITE(11,*) "# Cas B"
DO i=1,n

WRITE(11,*) zb(i),phib(i)

END DO
WRITE(11,"(/)") 

! b)

WRITE(11,*) "# Energia cas A"

WRITE(11,*) "Energia al principi(t=0): ",(lambda/2.d0)*za(1)**2 + omega*dsqrt(1-za(i)**2)*dcos(phia(1))&
,"Energia al final(t=10): ",(lambda/2.d0)*za(n)**2 + omega*dsqrt(1-za(i)**2)*dcos(phia(n))


WRITE(11,"(/)") 

WRITE(11,*) "# Energia cas B"

WRITE(11,*) "Energia al principi(t=0): ",(lambda/2.d0)*zb(1)**2 + omega*dsqrt(1-zb(i)**2)*dcos(phib(1))&
,"Energia al final(t=10): ",(lambda/2.d0)*zb(n)**2 + omega*dsqrt(1-zb(i)**2)*dcos(phib(n))


! 2)
! a)
X1 = 0.d0
x2 = 4.d0*datan(1.d0)
k = 13
CALL trapezoidalrule(x1,x2,k,f,resultat_tancat)

WRITE(11,"(/)") 

WRITE(11,*) "# Integral metode tancat"

WRITE(11,"(1e17.10)") resultat_tancat

! b)

WRITE(11,"(/)") 

WRITE(11,*) "# Convergencia de MonteCarlo amb mostreig d'importancia"

ndades = 100000
allocate(numeros(ndades))
xlow = 0.d0
xhigh = x2
cotasup = 1.5d0


CALL accepta(ndades,numeros,xlow,xhigh,cotasup,g)

DO i = 10000,100000,10000

CALL montecarlo(numeros,i,h,resultat,error)
WRITE(11,*) i,resultat,error
END DO

CLOSE(11)

END PROGRAM Exam

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

DOUBLE PRECISION :: x,lambda,omega
DOUBLE PRECISION ,DIMENSION(nequs)::funcin,dyoutput
INTEGER :: nequs
COMMON/CONSTANTS/lambda,omega

dyoutput(1) = -omega*dsqrt(1-funcin(1)**2)*dsin(funcin(2))
dyoutput(2) =  lambda*funcin(1) + omega*(funcin(1)/dsqrt(1-funcin(1)**2))*dcos(funcin(2))


END

SUBROUTINE integralRK4(x0,xf, N,nequs,funcin, z, phi)
IMPLICIT NONE
DOUBLE PRECISION :: x0,xf,dx,x
DOUBLE PRECISION,DIMENSION(nequs) :: funcin,phiRK
DOUBLE PRECISION,DIMENSION(N) ::  phi,z
INTEGER :: N, nequs, i
EXTERNAL EDO

dx = (xf-x0)/dble(N) 
DO i =1,n
x = x0 + dx*i
phi(i) = funcin(2)
z(i) = funcin(1)
call RungeKutta4(x,dx,funcin,phiRK,nequs,EDO)
funcin = phiRK
END DO

END

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

DOUBLE PRECISION FUNCTION f(x)

IMPLICIT NONE

DOUBLE PRECISION :: x,x0,gamma

gamma = 0.5d0
x0 = 1.d0

f = dsin(x)**2/((x-x0)**2 + gamma**2)

END FUNCTION

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

DOUBLE PRECISION FUNCTION g(x)

IMPLICIT NONE

DOUBLE PRECISION :: x


g = dsin(x)**2/(2.d0*datan(1.d0))

END FUNCTION

DOUBLE PRECISION FUNCTION h(x)

IMPLICIT NONE

DOUBLE PRECISION :: x,x0,gamma

gamma = 0.5d0
x0 = 1.d0

h = (2.d0*datan(1.d0))/((x-x0)**2 + gamma**2)

END FUNCTION