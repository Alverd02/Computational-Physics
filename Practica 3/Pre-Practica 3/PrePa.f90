PROGRAM P3

IMPLICIT NONE




END PROGRAM P3

SUBROUTINE newtonraphson(fun,xini,preci,nitera,valorarrel)

DOUBLE PRECISION :: xini,preci,valorarrel,x1,valordf,valorf
INTEGER :: nitera,
EXTERNAL fun

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

DOUBLE PRECISION :: A,B,preci,valorarrel,valorf_A,valorf_B
INTEGER :: nitera
EXTERNAL fun
CALL fun(A,valorf_A,valordf_A)
CALL fun(B,valorf_B,valordf_B)
IF ((valorf_A*valorf_B).GT.0.) THEN
WRITE(*,*) "No es pot aplicar el Teorema de Bolzano"
RETURN
END IF



RETURN 
END

SUBROUTINE fun(x,valorf,valordf)

DOUBLE PRECISION :: x,valorf,valordf

valorf = 
valordf = 

RETURN
END