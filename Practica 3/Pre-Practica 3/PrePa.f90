PROGRAM P3

IMPLICIT NONE

COMMON/CONSTANTS/pi 
DOUBLE PRECISION :: pi,valorarrel_1,fun,A,B,preci,valorarrel_2,valorarrel_3
INTEGER :: nitera_1,nitera_2,nitera_3
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

END PROGRAM P3
! Apartat 1

SUBROUTINE newtonraphson(fun,xini,preci,nitera,valorarrel)
IMPLICIT NONE
DOUBLE PRECISION :: xini,preci,valorarrel,x1,valordf,valorf
INTEGER :: nitera,i


DO i=1,999999

CALL fun(xini,valorf,valordf)
x1 = xini - valorf/valordf
IF (ABS(valorf/valordf).LE.preci) THEN
valorarrel = x1
nitera = i
EXIT
END IF
END DO

xini = x1

RETURN 
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
DOUBLE PRECISION :: x,valorf,valordf
COMMON/CONSTANTS/pi 
DOUBLE PRECISION :: pi

valorf = (-57./160*pi + (57./80 + 17./20*pi)*x - (17./10 + pi/2.)*x**2 + x**3)*sinh(x)
!valordf = ((57./80 + 17./20*pi) - 2*(17./10 + pi/2.)*x + 3*x**2)*sinh(x) + (-57./160*pi + (-57./160*pi + (57./80 + 17./20*pi)*x - (17./10 + pi/2.)*x**2 + x**3)*cosh(x)

END