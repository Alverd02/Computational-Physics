PROGRAM P6
! Al executar el programa se'm queda penjat
IMPLICIT NONE

INTEGER :: ISEED,i,ndades
DOUBLE PRECISION :: a,b,resultat_u,resultat_d,error_u,error_d,xlow,xhigh,cotasup,pi,error_2,I2,L,I3,error_3,u,d,p,g,f
DOUBLE PRECISION,DIMENSION(1000000) :: posis
COMMON/CONSTANTS/pi,L
EXTERNAL u
EXTERNAL d
EXTERNAL p
EXTERNAL g
EXTERNAL f

ISEED=20470586
CALL SRAND(ISEED)

pi = 4.d0*atan(1.d0)

OPEN(11,file="P6-23-24-res.dat")

! 1)
! a)
WRITE(11,*) "#0 N,n_u,error_u,n_d,error_d"

a = 0.d0
b = 1.d0

DO i=150,45000,150

CALL montecarlocru(a,b,i,u,resultat_u,error_u)
CALL montecarlocru(a,b,i,d,resultat_d,error_d)

WRITE(11,*) i,resultat_u,error_u,resultat_d,error_d

END DO

WRITE(11,"(/)")

! b)

L = pi
ndades = 1000000
xhigh = L 
xlow = -L
cotasup = 0.5

CALL accepta(ndades,posis,xlow,xhigh,cotasup,p)

! c)

WRITE(11,*) "#1 N,I2,ERROR_2"

DO i = 10000,1000000,10000

CALL montecarlo(posis,i,g,I2,error_2)

WRITE(11,*) i,I2,error_2

END DO

WRITE(11,*) "#2 N,I3,error_3"

xhigh = 1.d0
xlow = -1.d0
cotasup = 1.5

DO i=10000,300000,10000

CALL intaccepta(i,I3,error_3,xlow,xhigh,cotasup,f)
WRITE(11,*) i,2*I3,2*error_3
END DO

CLOSE(11)

END PROGRAM P6

DOUBLE PRECISION FUNCTION u(x)

IMPLICIT NONE

DOUBLE PRECISION :: x,u

u = (5.109*x**(0.8002)*(1-x)**3)/x

RETURN
END

DOUBLE PRECISION FUNCTION d(x)

IMPLICIT NONE

DOUBLE PRECISION :: x,d

d = (3.058*x**(0.803)*(1-x)**4)/x

RETURN
END

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

SUBROUTINE accepta(ndades,numeros,xlow,xhigh,cotasup,funcio)

IMPLICIT NONE

DOUBLE PRECISION :: xlow,xhigh,cotasup,funcio,p,valor_mitja,variancia,desv_est
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

DOUBLE PRECISION FUNCTION p(x)
IMPLICIT NONE
DOUBLE PRECISION :: x,p,pi,l
COMMON/CONSTANTS/pi,L

P = (1.d0/L)*(dsin((pi*(x-L))/(2*L)))**2

RETURN
END

SUBROUTINE montecarlo(valors,n,func,resultat,error)

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

DOUBLE PRECISION FUNCTION g(x)
 
IMPLICIT NONE 

DOUBLE PRECISION :: x,g,pi,l
COMMON/CONSTANTS/pi,L

g = (dsin((8*pi*(x-L))/(2*L)))**2

RETURN
END

SUBROUTINE intaccepta(ndades,resultat,error,xlow,xhigh,cotasup,funcio)

IMPLICIT NONE

DOUBLE PRECISION :: xlow,xhigh,cotasup,funcio,p,valor_mitja,variancia,desv_est,resultat,error
DOUBLE PRECISION, dimension(ndades) :: numeros
INTEGER :: ndades,n,i,n_dentro

n_dentro = 0

DO i=1,n


numeros(n) = (xhigh-xlow)*RAND() + xlow
p = cotasup*RAND()

IF (funcio(numeros(n)).ge.p) then

n_dentro = n_dentro + 1

END IF
END DO

resultat = (cotasup*(xhigh-xlow)*n_dentro)/n
error = ((cotasup*(xhigh-xlow))/dsqrt(dble(n)))*dsqrt((n_dentro/dble(n))*(1-n_dentro/n))
RETURN
END

DOUBLE PRECISION FUNCTION f(x)

IMPLICIT NONE

DOUBLE PRECISION :: x,f

f = dsqrt(1-x**2)

RETURN
END