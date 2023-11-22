       PROGRAM PREPRACTICA6
       IMPLICIT NONE
       EXTERNAL FUN1,FUN2,RHO,FUN3,FUN4,FUN5,FUN6

       INTEGER N,J
       DOUBLE PRECISION A,B,M,I1,I2,DESVEST1,DESVEST2,
     1 E,PI,REAL1,REAL2,ERROR1,ERROR2,FUN1,FUN2,RHO,MU,SIGMA,
     1 FUN3,FUN4,FUN5,DESVEST3,DESVEST4,DESVEST5,I3,I4,I5,
     1 XALRHO(1050000),XALGAUSS(1050000),I6,FUN6,DESVEST6
       COMMON/CONSTANTS/PI,E

       PI=4*DATAN(1.D0)
       E=DEXP(1.D0)


C 1)
C      A)

       OPEN(14,FILE="P6_1920_res.dat")
       WRITE(14,16) "N", "I1", "SIGMA1", "ERROR1", "I2", "SIGMA2", 
     1 "ERROR2"
16     FORMAT(4X,7(A,9X))
C valors reals de les integrals:
       REAL1=E*DSQRT((E**2)+(PI**2))+
     1 (PI**2)*DLOG((E/PI)+DSQRT(1+(E/PI)**2))
       REAL2=(1061.D0*PI/288)-(5.D0*PI**3)/12
C estimació per diferents N, amb la corresponent desviació i error
       DO N=2500,150000,2500
        CALL MONTECRU(N,-E,E,FUN1,I1,DESVEST1)
        ERROR1=DABS(REAL1-I1)
        CALL MONTECRU(N,-PI,PI,FUN2,I2,DESVEST2)
        ERROR2=DABS(REAL2-I2)
        WRITE(14,15) N, I1, DESVEST1, ERROR1, I2, DESVEST2, ERROR2
       ENDDO
15     FORMAT(I6,2X,E12.6,2X,2(E12.5,2X),E12.6,2X,2(E12.5,2X))

       CALL SYSTEM("gnuplot fig1.gnu")



C      B) i C)

       N=1050000
       M=1.3D0
       MU=0.D0
       SIGMA=1.D0
C guardo els nombres generats en un vector
       DO J=1,N,1
        CALL ACCEPTREGUIG(XALRHO(J),-PI,PI,M,RHO)
        CALL GAUSSIANA(XALGAUSS(J),MU,SIGMA)
       ENDDO
       
C      D)

C       WRITE(14,*)
C       WRITE(14,*)
C       WRITE(14,16) "N","I3","SIGMA3","I4","  SIGMA4",
C     1 "   I5","   SIGMA5"
C       DO J=5000,1050000,5000
C        CALL MONTECARLO(J,XALRHO,FUN3,I3,DESVEST3)
C        CALL MONTECARLO(J,XALRHO,FUN4,I4,DESVEST4)
C        CALL MONTECARLO(J,XALGAUSS,FUN5,I5,DESVEST5)
C        WRITE(14,19) J, I3, DESVEST3, I4, DESVEST4, I5,DESVEST5
C       ENDDO
C19     FORMAT(I10,2X,3(E12.6,2X,E12.5,2X))

c manera més eficient:
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "N   "," I3    ","  SIGMA3"
       CALL MONTEFIT(1050000,XALRHO,FUN3,I3,DESVEST3)
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "N   "," I4    ","  SIGMA4"
       CALL MONTEFIT(1050000,XALRHO,FUN4,I4,DESVEST4)
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "N   "," I5    ","  SIGMA5"
       CALL MONTEFIT(1050000,XALGAUSS,FUN5,I5,DESVEST5)



C 2)
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "N   "," I6    ","  SIGMA6"
       CALL MONTEMULTI(1050000,XALGAUSS,FUN6,I6,DESVEST6)

       CLOSE(14)


       END PROGRAM



C DEFINICIÓ SUBRUTINES


C Montecarlo cru
       SUBROUTINE MONTECRU(N,A,B,FUNCI,I,DESVEST)
       IMPLICIT NONE
       INTEGER N,J
       DOUBLE PRECISION FUNCI,I,DESVEST,X,H,A,B,HQUAD,IQUAD
       I=0.D0
       IQUAD=0.D0
       DO J=1,N,1
C genera un nombre aleatori entre 0 i 1
        X=RAND()
C en calcula el valor de la funció amb el canvi de variable i el seu quadrat per després calcular la desviació estàndard
        H=(B-A)*FUNCI((B-A)*X+A)
        HQUAD=H**2
        I=I+H
        IQUAD=IQUAD+HQUAD
       ENDDO
C valor final de la integral:
       I=I/N
       IQUAD=IQUAD/N
C desviació estàndard:
       DESVEST=DSQRT((IQUAD-I**2)/N)
       RETURN
       END


C mètode d'acceptació-rebuig
       SUBROUTINE ACCEPTREGUIG(XAL,A,B,M,FUN)
       IMPLICIT NONE
       DOUBLE PRECISION A,B,XAL,M,X,P,X1,X2,FUN
C genera x1 i x2 aleatoris dins de U(0,1)
21     X1=RAND()
       X2=RAND()
C canvi de variable perquè x pertanyi a U(a,b) i p a U(0,M)
       X=(B-A)*X1+A
       P=M*X2
C condició: fun(x)>=p
       IF (FUN(X).LT.P) THEN
        GO TO 21
       ENDIF
       XAL=X
       RETURN
       END


C mètode de Box-Müller
       SUBROUTINE GAUSSIANA(XAL,MU,SIGMA)
       IMPLICIT NONE
       DOUBLE PRECISION XAL,MU,SIGMA,PI,R,PHI
       PI=4*DATAN(1.D0)
C genera els nombres r i phi a partir de U(0,1) amb un canvi de variable
       R=DSQRT(-2.D0*LOG(RAND()))
       PHI=2*PI*RAND()
C nombre aleatori segons N(x,mu,sigma):
       XAL=MU+SIGMA*(R*DCOS(PHI))
       RETURN
       END
       

C Montecarlo d'importància
       SUBROUTINE MONTECARLO(N,XAL,FUNCI,I,DESVEST)
       IMPLICIT NONE
       INTEGER N,J
       DOUBLE PRECISION XAL(N),X,FUNCI,I,DESVEST,IQUAD
       I=0.D0
       IQUAD=0.D0
       DO J=1,N,1
C en calcula el valor de la funció i el seu quadrat per després calcular la desviació estàndard
        I=I+FUNCI(XAL(J))
        IQUAD=IQUAD+(FUNCI(XAL(J)))**2
       ENDDO
C valor final de la integral:
       I=I/N
       IQUAD=IQUAD/N
C desviació estàndard:
       DESVEST=DSQRT((IQUAD-I**2)/N)
       RETURN
       END


C Montecarlo d'importància que guarda resultats al fitxer cada N=5000
c Calcula la integral I de la funció FUNCI i la seva desviació estàndard DESVEST
C amb els nombres aleatoris del vector XAL(N).
c Només fa el càlcul 1 cop fins la N màxima i guarda resultats parcials.
       SUBROUTINE MONTEFIT(N,XAL,FUNCI,I,DESVEST)
       IMPLICIT NONE
       INTEGER N,J
       DOUBLE PRECISION XAL(N),X,FUNCI,I,DESVEST,IQUAD,
     1 IPARCIAL,IQUADPAR
       I=0.D0
       IQUAD=0.D0
       DO J=1,N,1
C en calcula la integral acumulada i el seu quadrat
        I=I+FUNCI(XAL(J))
        IQUAD=IQUAD+(FUNCI(XAL(J)))**2
        IF (MOD(J,5000).EQ.0) THEN
          IPARCIAL=I/J
          IQUADPAR=IQUAD/J
C desviació estàndard:
          DESVEST=DSQRT((IQUADPAR-IPARCIAL**2)/J)
         WRITE(14,23) J, IPARCIAL, DESVEST
        ENDIF
       ENDDO
23     FORMAT(I10,2X,E12.6,2X,E12.5)
       RETURN
       END



C Montecarlo multidimensional
       SUBROUTINE MONTEMULTI(N,XAL,FUNCI,I,DESVEST)
       IMPLICIT NONE
       INTEGER N,J
       DOUBLE PRECISION XAL(N),FUNCI,I,DESVEST,IQUAD,
     1 IPARCIAL,IQUADPAR
       I=0.D0
       IQUAD=0.D0
       DO J=5,N,5
C en calcula la integral acumulada i el seu quadrat
        I=I+FUNCI(XAL(J),XAL(J-1),XAL(J-2),XAL(J-3),XAL(J-4))
        IQUAD=IQUAD+(FUNCI(XAL(J),XAL(J+1),XAL(J+2),XAL(J+3),
     1 XAL(J+4)))**2
C guarda resultats parcials cada 1500 sumands
        IF (MOD(J,1500*5).EQ.0) THEN
          IPARCIAL=I*5.D0/J
          IQUADPAR=IQUAD*5.D0/J
C desviació estàndard:
         DESVEST=DSQRT((IQUADPAR-IPARCIAL**2)*5.D0/J)
         WRITE(14,23) J/5, IPARCIAL, DESVEST
        ENDIF
       ENDDO
23     FORMAT(I10,2X,E12.6,2X,E12.5)
       RETURN
       END



C DEFINICIÓ DE LES FUNCIONS

C funció a integrar 1
       DOUBLE PRECISION FUNCTION FUN1(X)
       IMPLICIT NONE
       DOUBLE PRECISION PI,X,E
       COMMON/CONSTANTS/PI,E
       FUN1=DSQRT((PI**2)+(X**2))
       RETURN
       END

C funció a integrar 2
       DOUBLE PRECISION FUNCTION FUN2(X)
       IMPLICIT NONE
       DOUBLE PRECISION PI,X,E
       COMMON/CONSTANTS/PI,E
       FUN2=(X+3*(X**2)*DSIN(X)-(X**3))*((DCOS(X))**2)*DSIN(X)
       RETURN
       END

C distribució de probabilitat
       DOUBLE PRECISION FUNCTION RHO(X)
       IMPLICIT NONE
       DOUBLE PRECISION PI,E,X
       COMMON/CONSTANTS/PI,E
       RHO=((5.D0/4)*DEXP(-DABS(X))*(DSIN(X))**2)/
     1 (1-DEXP(-PI))
       RETURN
       END


C funció a integrar 3
       DOUBLE PRECISION FUNCTION FUN3(X)
       IMPLICIT NONE
       DOUBLE PRECISION X,PI,E
       COMMON/CONSTANTS/PI,E
C funció a integrar dividida per la distribució de probabilitat rho(x)
       FUN3=4.D0*(1-DEXP(-PI))*(X**2)/5
       RETURN
       END

C funció a integrar 4
       DOUBLE PRECISION FUNCTION FUN4(X)
       IMPLICIT NONE
       DOUBLE PRECISION X,PI,E
       COMMON/CONSTANTS/PI,E
C funció a integrar dividida per la distribució de probabilitat rho(x)
       FUN4=4.D0*(1-DEXP(-PI))*(1+X**2)*
     1 DEXP(DABS(X)-(X**2)/2)/(5*(DTAN(X))**2)
       RETURN
       END

C funció a integrar 5
       DOUBLE PRECISION FUNCTION FUN5(X)
       IMPLICIT NONE
       DOUBLE PRECISION X,PI,E
       COMMON/CONSTANTS/PI,E
C funció a integrar dividida per la distribució de probabilitat gaussiana
       FUN5=DSQRT(2*PI)*DEXP(X**2/2)*(DSIN(X))**4*X**2
       RETURN
       END

C funció a integrar 6
       DOUBLE PRECISION FUNCTION FUN6(X1,X2,X3,X4,X5)
       IMPLICIT NONE
       DOUBLE PRECISION X1,X2,X3,X4,X5,PI,E
       COMMON/CONSTANTS/PI,E
C funció a integrar dividida per la distribució de probabilitat gaussiana de les 5 variables
       FUN6=(2*PI)**(-5)*DEXP(-(X1**2+X2**2+X3**2+X4**2+X5**2)/2)*
     1 DEXP(X1*DCOS(X2+X3))*(X3**2*X4**2*X5**2+(DCOS(X3+X4))**2*
     1 DSIN(X5))
       RETURN
       END
