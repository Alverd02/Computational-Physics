PROGRAM P6

IMPLICIT NONE

INTEGER :: ISEED,i,ndades,j,K,n
DOUBLE PRECISION :: pi,e,resultat1,resultat2,error1,error2,a,b,cotasup,sigma,mu,I1,I2
DOUBLE PRECISION, DIMENSION(1100000) :: numeros,yres
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

I1 = pi**3/2.d0
I2 = 6*pi**3 - (905/144)*pi

OPEN(11,file="P6-23-24-res.dat")

WRITE(11,*) "# i,res1,sigma1,err1,res2,sigma2,err2"

DO i=2000,120000,2000
CALL montecarlocru(-pi,pi,i,func1,resultat1,error1)
CALL montecarlocru(-2*pi,2*pi,i,func2,resultat2,error2)

WRITE(11,*) i,resultat1,error1,abs(resultat1-I1),resultat2,error2,abs(resultat2-I2)

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


WRITE(11,"(/)")
CLOSE(11)





END PROGRAM P6

subroutine montecarlocru(a,b,n,func,resultat,error)

IMPLICIT NONE

INTEGER :: n,i
DOUBLE PRECISION :: h,a,b,func,resultat,error,x2
DOUBLE PRECISION, DIMENSION(n) :: t

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