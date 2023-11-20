PROGRAM P6

IMPLICIT NONE

INTEGER :: ISEED,i,ndades,k,j
DOUBLE PRECISION :: pi,e,resultat1,resultat2,error1,error2,a,b,cotasup,sigma,mu,resultat3,resultat4,resultat5,error5,error4,error3
DOUBLE PRECISION, DIMENSION(1100000) :: numeros,yres,t2
DOUBLE PRECISION, DIMENSION(120000) :: t
COMMON/CONSTANTS/pi,e
EXTERNAL func1
EXTERNAL func2
EXTERNAL p
EXTERNAL func3
EXTERNAL func4
EXTERNAL func5



ISEED=20470586
CALL SRAND(ISEED)

pi = 3.1415926535898 
e = 2.718281828459 

OPEN(11,file="P6-23-24-res.dat")

WRITE(11,*) "# a)"

DO i=2000,120000,2000

DO k = 1,i

t(k) = RAND()

END DO

CALL montecarlo(-pi,pi,i,func1,resultat1,error1,t)
CALL montecarlo(-2*pi,2*pi,i,func2,resultat2,error2,t)

WRITE(11,*) i,resultat1,error1,resultat2,error2

END DO

WRITE(11,"(/)")

ndades = 1100000
a = 0
b = pi
cotasup = ((10/3.d0)*e**(-pi/2.d0))/(1 + e**(-pi))
CALL accepta(ndades,numeros,a,b,cotasup,p)

sigma = 1/dsqrt(2.d0)
mu = 0.d0

CALL boxmuller(ndades,sigma,mu,yres)

WRITE(11,*) "# d)"

DO i = 5000,1100000,5000

DO k = 1,i

t2(k) = RAND()

END DO

CALL montecarlo(0,pi,i,func3,resultat3,error3,t2)
CALL montecarlo(0,pi,i,func4,resultat4,error4,t2)
!CALL montecarlo(-9999999,9999999,i,func5,resultat5,error5,yres)
  
WRITE(11,*) i,resultat3,error3,resultat4,error4,resultat5,error5

END DO
WRITE(11,"(/)")
CLOSE(11)





END PROGRAM P6

subroutine montecarlo(a,b,n,func,resultat,error,t)

IMPLICIT NONE

INTEGER :: n,i
DOUBLE PRECISION :: h,a,b,func,resultat,error,x2
DOUBLE PRECISION, DIMENSION(n) :: t

resultat = 0.d0
x2 = 0.d0

DO i=1,n 

h = (b-a)*func((b-a)*t(i)+a)
resultat = resultat + h
x2 = x2 + h**2

END DO

error = (1/dsqrt(dble(n)))*dsqrt((x2/dble(n)) - (resultat/dble(n))**2)
resultat = resultat/dble(n)

RETURN
END

DOUBLE PRECISION FUNCTION func1(x)
IMPLICIT NONE
COMMON/CONSTANTS/pi
DOUBLE PRECISION :: pi,x,y

y = dsqrt(pi**2-x**2)

RETURN
END

DOUBLE PRECISION FUNCTION func2(x)
IMPLICIT NONE
DOUBLE PRECISION :: x,y

y = (x**2*dsin(x)-x**3)*dcos(x)**2*dsin(x)

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

DOUBLE PRECISION :: x,y,pi,e
COMMON/CONSTANTS/pi,e

y = ((10/3.d0)*e**(-x)*(sin(x))**3)/(1+e**(-pi))

RETURN
END

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

DOUBLE PRECISION FUNCTION func3(x)

IMPLICIT NONE

DOUBLE PRECISION :: x,y,pi,e
COMMON/CONSTANTS/pi,e

y = e**(-abs(x))*x**2*(dsin(x))**2

RETURN
END

DOUBLE PRECISION FUNCTION func4(x)

IMPLICIT NONE

DOUBLE PRECISION :: x,y,pi,e
COMMON/CONSTANTS/pi,e

y = e**(-(x**2/2.d0))*(cos(x))**2*(pi + 4*x**2)

RETURN
END

DOUBLE PRECISION FUNCTION func5(x)

IMPLICIT NONE

DOUBLE PRECISION :: x,y,pi,e
COMMON/CONSTANTS/pi,e

y = e**(-x**2)*(dsin(x))**2*x**2

RETURN
END