PROGRAM P3

IMPLICIT NONE

COMMON/CONSTANTS/pi 
DOUBLE PRECISION :: pi,valorarrel_1,A,B,preci,valorarrel_2,valorarrel_3,xini,valorarrel,E,valordf,valorf
DOUBLE PRECISION, DIMENSION(1:9) :: E_0
DOUBLE PRECISION, DIMENSION(1:25) :: E_25,funci_25,dfunci_25,dfunci_25_aprox
DOUBLE PRECISION, DIMENSION(1:230) :: E_230,funci_230,dfunci_230,dfunci_230_aprox
INTEGER :: nitera_1,nitera_2,nitera_3,nitera,i,ndates
EXTERNAL fun
pi = 3.1415926535898


! Apartat 2
A = 1.5d0
B = 2.0d0
preci = 1.d-12

CALL biseccio(fun,A,B,preci,nitera_1,valorarrel_1)
WRITE(*,"(A21,I3,A18,e20.12)") "Nombre d'iteracions 1: ",nitera_1," Valor de l'arrel 1: ",valorarrel_1

A = 0.5d0
B = 0.8d0

CALL biseccio(fun,A,B,preci,nitera_2,valorarrel_2)
WRITE(*,"(A21,I3,A18,e20.12)") "Nombre d'iteracions 2: ",nitera_2," Valor de l'arrel 2: ",valorarrel_2

A = 0.8d0
B = 1.0d0

CALL biseccio(fun,A,B,preci,nitera_3,valorarrel_3)
WRITE(*,"(A21,I3,A18,e20.12)") "Nombre d'iteracions 3: ",nitera_3," Valor de l'arrel 3: ",valorarrel_3

E_0 = [0.1,0.2,0.65,0.7,1.3,2.4,2.6,3.9,5.3]
OPEN(12,file="P3-23-24-res.dat")
DO i=1,size(E_0)

xini = E_0(i)

CALL newtonraphson(fun,xini,preci,nitera,valorarrel)

END DO
CLOSE(12)

! Apartat 4

DO i=0,25
E = 0.251327*i
E_25(i) = E
CALL fun(E,valorf,valordf)
funci_25(i) = valorf
dfunci_25(i) = valordf
END DO

DO i=0,230
E = 0.027318*i
E_230(i) = E
CALL fun(E,valorf,valordf)
funci_230(i) = valorf
dfunci_230(i) = valordf
END DO


OPEN(13,file =" P3-23-24-res3-n25.dat")
ndates = 25
CALL derivataula(ndates,E_25,funci_25,dfunci_25_aprox)

DO i=1,ndates

WRITE(13,"(4e20.12))") E_25(i),funci_25(i),dfunci_25_aprox(i),dfunci_25(i)
END DO

CLOSE(13)

OPEN(14,file =" P3-23-24-res3-n230.dat")
ndates = 230
CALL derivataula(ndates,E_230,funci_230,dfunci_230_aprox)

DO i=1,ndates

WRITE(14,"(4e20.12))") E_230(i),funci_230(i),dfunci_230_aprox(i),dfunci_230(i)
END DO
CLOSE(14)

END PROGRAM P3
! Apartat 1

SUBROUTINE newtonraphson(fun,xini,preci,nitera,valorarrel)
IMPLICIT NONE
DOUBLE PRECISION :: xini,preci,valorarrel,x1,valordf,valorf
INTEGER :: nitera,i

write(12, "(A23,f4.2)") "# Resultados para E0 = ", xini
DO i=1,999999

CALL fun(xini,valorf,valordf)
x1 = xini - valorf/valordf
IF (ABS(valorf/valordf).LE.preci) THEN
valorarrel = x1
nitera = i
WRITE(12,"(1I2,e20.12)") nitera,valorarrel
WRITE(12,"(/)")
RETURN 
STOP
END IF
valorarrel = x1
nitera = i

WRITE(12,"(1I2,e20.12)") nitera,valorarrel
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

SUBROUTINE fun(x,valorf,valordf)
IMPLICIT NONE
DOUBLE PRECISION :: x,valorf,valordf,valordf_1,valordf_2
COMMON/CONSTANTS/pi 
DOUBLE PRECISION :: pi

valorf = (-57./160*pi + (57./80 + 17./20*pi)*x - (17./10 + pi/2.)*x**2 + x**3)*sinh(x)
valordf_1 = ((57./80 + 17./20*pi) - 2*(17./10 + pi/2.)*x + 3*x**2)*sinh(x) 
valordf_2 = (-57./160*pi + (57./80 + 17./20*pi)*x - (17./10 + pi/2.)*x**2 + x**3)*cosh(x)
valordf = valordf_1 + valordf_2
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