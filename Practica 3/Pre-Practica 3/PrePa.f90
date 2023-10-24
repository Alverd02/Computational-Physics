PROGRAM P3

IMPLICIT NONE

COMMON/CONSTANTS/pi 
DOUBLE PRECISION :: pi,valorarrel_1,fun
INTEGER :: nitera_1
EXTERNAL fun
pi = 3.1415926535898

! Apartat 2

CALL biseccio(fun,1.5d0,2.0d0,1.d-12,nitera_1,valorarrel_1)
WRITE(*,"(e20.12)") nitera_1,valorarrel_1
END PROGRAM P3

! Apartat 1

SUBROUTINE newtonraphson(fun,xini,preci,nitera,valorarrel)

DOUBLE PRECISION :: xini,preci,valorarrel,x1,valordf,valorf
INTEGER :: nitera


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

DOUBLE PRECISION :: A,B,preci,valorarrel,valorf_A,valorf_B,d,valordf_A,valordf_B,valorf_d,valordf_d
INTEGER :: i,nitera

nitera=INT(LOG((B-A)/preci)/LOG(2.0d0))+1
DO i=1,nitera
d = (A+B)/2
CALL fun(A,valorf_A,valordf_A)
CALL fun(B,valorf_B,valordf_B)
CALL fun(d,valorf_d,valordf_d)


IF ((valorf_A*valorf_d).LT.0.D0) THEN

B = d

ELSE

A = d

END IF
END DO

valorarrel = d
nitera=INT(LOG((B-A)/preci)/LOG(2.0d0))+1
RETURN 
END

SUBROUTINE fun(x,valorf,valordf)

DOUBLE PRECISION :: x,valorf,valordf
COMMON/CONSTANTS/pi 
DOUBLE PRECISION :: pi

valorf = (-57./160*pi + (57./80 + 17./20*pi)*x - (17./10 + pi/2.)*x**2 + x**3)*sinh(x)
!valordf = ((57./80 + 17./20*pi) - 2*(17./10 + pi/2.)*x + 3*x**2)*sinh(x) + (-57./160*pi + (-57./160*pi + (57./80 + 17./20*pi)*x - (17./10 + pi/2.)*x**2 + x**3)*cosh(x)

END