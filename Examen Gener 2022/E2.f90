PROGRAM E2

IMPLICIT NONE
DOUBLE PRECISION :: sigma,A,M,C,T_ext,P,dt,tf,preci,T_sat,f,df
INTEGER :: k,nitera
common/CONSTANTS/sigma,A,M,C,T_ext,P
EXTERNAL dt
EXTERNAL f
EXTERNAL df

sigma = 5.67d-8
A = 2.04387d0
M = 70.d0
    C = 4184.d0
T_ext = 2.3d0
P = 0.d0

OPEN(11,file="apartat2a.dat")

k = 10
CALL trapezoidalrule(310.d0,304.d0,k,dt,tf)
WRITE(11,*) "# Temps de refredament"
WRITE(11,*) tf
WRITE(11,"(/)")

P = 100.D0
CALL trapezoidalrule(310.d0,304.d0,k,dt,tf)
WRITE(11,*) "# Temps de refredament"
WRITE(11,*) tf
WRITE(11,"(/)")


CLOSE(11)

OPEN(12,file="apartat2b.dat")
WRITE(12,*) "#No fonts"
P = 0
CALL trapezoidalrule(310.d0,304.d0,k,dt,tf)
WRITE(12,"(/)")
WRITE(12,*) "#Fonts"
P = 100.D0
CALL trapezoidalrule(310.d0,304.d0,k,dt,tf)
CLOSE(12)


preci = 1D-8
P = 100D0
call newtonraphson(f,df,200.d0,preci,nitera,T_sat)
write(*,*) "Temperatura de Saturación para el caso con fuente (K)"
write(*,*) T_sat


END PROGRAM E2

DOUBLE PRECISION FUNCTION dt(T)

IMPLICIT NONE
DOUBLE PRECISION :: sigma,A,M,C,T_ext,P,dt,T
common/CONSTANTS/sigma,A,M,C,T_ext,P



dt = (M*C)/(sigma*A*(T_ext**4-T**4) + P)
RETURN
END FUNCTION

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
x = x1+i*h
valor =  func(x)
resultat2 = valor + resultat2
resultat = resultat1 + h*resultat2
WRITE(12,*) x,resultat
END DO

resultat = resultat1 + h*resultat2
RETURN

END

SUBROUTINE newtonraphson(f,df,xini,preci,nitera,valorarrel)
IMPLICIT NONE
DOUBLE PRECISION :: xini,preci,valorarrel,x1,valordf,valorf,f,df
INTEGER :: nitera,i


DO i=1,999999

valorf = f(xini)  
valordf = df(xini)  
x1 = xini - valorf/valordf


IF (ABS(valorf/valordf).LE.preci) THEN
valorarrel = x1
nitera = i
RETURN 
STOP
END IF

nitera = i

xini = x1

END DO

END

DOUBLE PRECISION FUNCTION f(T)
IMPLICIT NONE
DOUBLE PRECISION :: sigma,A,M,C,T_ext,P,T,f
common/CONSTANTS/sigma,A,M,C,T_ext,P

f = (sigma*A*(T_ext**4-T**4) + P)/(M*C)

RETURN
END

DOUBLE PRECISION FUNCTION df(T)
IMPLICIT NONE
DOUBLE PRECISION :: sigma,A,M,C,T_ext,P,T,df
common/CONSTANTS/sigma,A,M,C,T_ext,P

df = -(sigma*A*4.d0*T**3)/(M*C)

RETURN
END