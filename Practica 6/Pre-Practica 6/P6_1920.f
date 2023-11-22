C FERRAN SEBARROJA TORRA, M1
C Física computacional, 2019/20
C Pràctica 6: Nombres aleatoris 2


       PROGRAM PREPRACTICA6
       IMPLICIT NONE
       EXTERNAL U,D,RHO,G,PHI

       INTEGER ISEED,N,J
       DOUBLE PRECISION U,D,X,NU,ND,SIGMANU,SIGMAND,RHO,PI,
     1 XALRHO(1000000),I2,SIGMA2,L,G,PHI,I3,SIGMA3

       PI=4*DATAN(1.D0)

       ISEED=20034486
       CALL SRAND(ISEED)

C 1)
C      a)
       OPEN(14,FILE="P6_1920_res.dat")
       WRITE(14,*) "N,   Nu,   SIGMAu,   Nd,    SIGMAd"
       DO N=150,45000,150
        CALL MONTECRU(N,0.D0,1.D0,U,NU,SIGMANU)
        CALL MONTECRU(N,0.D0,1.D0,D,ND,SIGMAND)
        WRITE(14,15) N,NU,SIGMANU,ND,SIGMAND
       ENDDO
15     FORMAT(I6,2X,2(E12.6,2X,E12.5,2X))

C      b)
C crea el vector de nombres aleatoris distribuïts segons RHO
       L=PI
       DO J=1,1000000,1
        CALL ACCEPTREGUIG(XALRHO(J),-L,L,1.D0/L,RHO)
       ENDDO

C      c)
c Calcula la integral amb els nombres aleatoris de RHO, guardant
C resultats parcials en un fitxer cada 10000 sumands.
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "N   "," I2    ","  SIGMA2"
       CALL MONTEFIT(1000000,XALRHO,G,I2,SIGMA2)


C 2)
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "N   "," I3    ","  SIGMA3"
       CALL MONTEMULTI(750000,XALRHO,PHI,I3,SIGMA3)



       CLOSE(14)

       END PROGRAM




C DEFINICIÓ SUBRUTINES--------------------------------------------------------


C Montecarlo cru
       SUBROUTINE MONTECRU(N,A,B,FUNCI,I,DESVEST)
C calcula la integral I i la desviació estandard DESVEST de la funció FUNCI
C entre A i B amb N nombres aleatoris.
       IMPLICIT NONE
       INTEGER N,J
       DOUBLE PRECISION FUNCI,I,DESVEST,X,H,A,B,HQUAD,IQUAD
       I=0.D0
       IQUAD=0.D0
       DO J=1,N,1
C genera un nombre aleatori entre 0 i 1
        X=RAND()
C en calcula el valor de la funció amb el canvi de variable i el seu quadrat 
c per després calcular la desviació estàndard
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
C retorna un nombre aleatori XAL entre A i B distribuït segons 
C la funció FUN, amb la cota màxima M
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


C Montecarlo d'importància que guarda resultats al fitxer cada N=10000
       SUBROUTINE MONTEFIT(N,XAL,FUNCI,I,DESVEST)
c Calcula la integral I de la funció FUNCI i la seva desviació estàndard DESVEST
C amb els nombres aleatoris del vector XAL(N). Només fa el càlcul 1 cop fins la 
C N màxima i guarda en el fitxer ober a 14 resultats parcials.
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
C guarda resultats parcials per múltiples de 10000
        IF (MOD(J,10000).EQ.0) THEN
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
C Calcula la integral I de la funció FUNCI de 3 variables amb els nombres 
c aleatoris del vector XAL(N), i la seva sigma DESVEST. Escriu el el fitxer 
C obert a 14 els resultats parcials cada N=10000.
       IMPLICIT NONE
       INTEGER N,J
       DOUBLE PRECISION XAL(N),FUNCI,I,DESVEST,IQUAD,
     1 IPARCIAL,IQUADPAR
       I=0.D0
       IQUAD=0.D0
       DO J=3,N,3
C en calcula la integral acumulada i el seu quadrat
        I=I+FUNCI(XAL(J),XAL(J-1),XAL(J-2))
        IQUAD=IQUAD+(FUNCI(XAL(J),XAL(J-1),XAL(J-2)))**2
C guarda resultats parcials cada 1000 sumands 
C (#sumands=J/3 amb J=# de nombres aleatoris utilitzats)
        IF (MOD(J,1000*3).EQ.0) THEN
          IPARCIAL=I*3.D0/J
          IQUADPAR=IQUAD*3.D0/J
C desviació estàndard:
         DESVEST=DSQRT((IQUADPAR-IPARCIAL**2)*3.D0/J)
         WRITE(14,23) J/3, IPARCIAL, DESVEST
        ENDIF
       ENDDO
23     FORMAT(I10,2X,E14.6,2X,E12.5)
       RETURN
       END







C DEFINICIÓ DE LES FUNCIONS

C densitat "up"
       DOUBLE PRECISION FUNCTION U(X)
       IMPLICIT NONE
       DOUBLE PRECISION X
       U=5.109*X**(0.8002-1)*(1-X)**3
       RETURN
       END


C densitat "down"
       DOUBLE PRECISION FUNCTION D(X)
       IMPLICIT NONE
       DOUBLE PRECISION X
       D=3.058*X**(0.803-1)*(1-X)**4
       RETURN
       END

C posició àtom ultrafred
       DOUBLE PRECISION FUNCTION RHO(X)
       IMPLICIT NONE
       DOUBLE PRECISION X,L,PI
       PI=4*DATAN(1.D0)
       L=PI
       RHO=L**(-1)*(DSIN(PI*(X-L)/(2*L)))**2
       RETURN
       END

C funció a integrar per I2
       DOUBLE PRECISION FUNCTION G(X)
       IMPLICIT NONE
       DOUBLE PRECISION X,PI,L
       PI=4*DATAN(1.D0)
       L=PI
       G=(DSIN(8*PI*(X-L)/(2*L)))**2
       RETURN
       END

C Mòdul al quadrat de la funció d'ona (3 variables)
C dividida per la distribució de les 3 variables
       DOUBLE PRECISION FUNCTION PHI(X1,X2,X3)
       IMPLICIT NONE
       INTEGER I,J,K
       DOUBLE PRECISION X1,X2,X3,PI,L,X(3),RHO(3)
       PI=4*DATAN(1.D0)
       L=PI
       X=[X1,X2,X3]
C productoris per I,J,K
       DO I=1,3,1
        DO K=1,3,1
         IF ((K.GT.1).AND.(K.LT.3)) THEN
          DO J=1,K-1,1
           PHI=(DSIN(PI*(X(I)-L)/(2*L))*
     1     (DCOS(PI*(X(J)-L)/(2*L))-
     1     DCOS(PI*(X(K)-L)/(2*L))))**2
          ENDDO
         ENDIF
        ENDDO
C aprofito el bucle de I entre 1 i 3 per calcular la distribució
C de cada variable i dividir-la a la funció PHI que s'integrarà.
       RHO(I)=L**(-1)*(DSIN(PI*(X(I)-L)/(2*L)))**2
       PHI=PHI/RHO(I)
       ENDDO
       RETURN
       END
