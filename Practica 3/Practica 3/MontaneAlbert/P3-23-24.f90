PROGRAM P3

IMPLICIT NONE

COMMON/CONSTANTS/pi 
DOUBLE PRECISION :: a,eps,E_i,D,pi,x,y,A_2,B,preci,valorarrel,t,xini,valordf,valorf
DOUBLE PRECISION, DIMENSION(1:100) :: E,posis_halley,dfunci_halley
INTEGER :: i,ndates,nitera
EXTERNAL fun_1
EXTERNAL fun_2
pi = 3.1415926535898

! Apartat 1

a = 17.857619
eps =  0.967990

DO i= 0,99

E(i+1) = 0.062831853*i

END DO

DO i=1,100

E_i = E(i)
CALL posis(E_i,a,eps,D,x,y)

posis_halley(i) = D

END DO

OPEN(11,file ="P3-23-24-res.dat")

ndates = 100
CALL derivataula(ndates,E,posis_halley,dfunci_halley)
WRITE(11,"(A23)") "# Resultados apartado 1"
DO i=1,ndates
WRITE(11,"(3e20.12))") E(i),posis_halley(i),dfunci_halley(i)
END DO
WRITE(11,"(/)")

! Apartat 2

A_2 = 0.1
B = 5.8
preci = 1.d-10

CALL biseccio(fun_1,A_2,B,preci,nitera,valorarrel)
CALL posis(valorarrel,a,eps,D,x,y)

WRITE(11,"(A23)") "# Resultados apartado 2"
WRITE(11,"(3e20.12))") valorarrel,D
WRITE(11,"(/)")

WRITE(11,"(A23)") "# Resultados apartado 3"

DO i=0,79
t = 0.94125*i
xini = pi/4.d0
preci = 1.d-12

CALL newtonraphson(fun_2,xini,preci,nitera,valorarrel)
CALL posis(valorarrel,a,eps,D,x,y)
WRITE(11,"(4e20.12))") t,valorarrel,x,y
END DO

CLOSE(11)

END PROGRAM P3

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

SUBROUTINE fun_1(E,valorf,eps)
IMPLICIT NONE
DOUBLE PRECISION :: E,valorf,eps,valorf_1,valorf_2

valorf_1 = sin(2*E)*(1-eps**2)
valorf_2 = -(cos(E)*(2-eps**2)-eps)*sin(E)
valorf = valorf_1 + valorf_2
RETURN
END

SUBROUTINE fun_2(E,valorf,eps,valordf,t)
IMPLICIT NONE
DOUBLE PRECISION :: E,valorf,eps,valordf,pi,T_H,t
COMMON/CONSTANTS/pi
T_H = 75.3
valorf = (2*pi*t)/(T_H)-eps*sin(E)
valordf = -eps*cos(E)

RETURN
END

! Apartat 3

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

SUBROUTINE posis(E,a,eps,D,x,y)
IMPLICIT NONE
DOUBLE PRECISION :: a,eps,D,y,x,E

x = a*(cos(E)-eps)
y = a*sqrt(1-eps**2)*sin(E)

D = sqrt(x**2+y**2)

RETURN
END