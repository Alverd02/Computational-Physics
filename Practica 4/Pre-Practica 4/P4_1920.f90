       PROGRAM PRE4
       IMPLICIT NONE
       EXTERNAL LONGITUD,DENSITAT,F3

       INTEGER K
       REAL XKR
       DOUBLE PRECISION XKDP,PI,RESULTATT,RESULTATS,L,H
       DOUBLE PRECISION EXACTEA,EXACTEM
       PARAMETER(PI=4.D0*DATAN(1.D0))


C 0)
C   a)
       XKR=0.
       XKDP=0.D0
       OPEN(14,FILE="dades.dat")
       WRITE(14,*) "Suma de reals i doble precisió"
       DO K=0,200000000,1
        XKR=XKR+0.02
        XKDP=XKDP+0.02D0
C condició que k sigui divisible per 100000 (el residu sigui 0)
        IF (MOD(K,100000).EQ.0) THEN
        WRITE(14,15) XKR, XKDP
        ENDIF
       ENDDO
15     FORMAT(F16.8,2X,F20.12)

C   b)
       XKR=0.
       XKDP=0.D0
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "Multiplicació de reals i doble precisió"
       DO K=0,2000,1
        XKR=2000.*K
        XKDP=2000.D0*K
        WRITE(14,15) XKR, XKDP
       ENDDO
       CLOSE(14)

       PRINT*, "Sumant es genera molt error, multiplicant no."
       PRINT*, "Si es multiplica tant és real com doble precisió."

C 1) Creació subrutines

C 2)
       OPEN(16,FILE="P4_1920_res1.dat")
C   a) Àrea sota una corba:
       CALL TRAPEZOIDS(-PI,PI,18,LONGITUD,RESULTATT)
       CALL SIMPSON(-PI,PI,18,LONGITUD,RESULTATS)
       WRITE(16,*) "Àrea calculada per trapezis i Simpson:"
       WRITE(16,*) RESULTATT,RESULTATS

C   b) Massa integrant la densitat:
       WRITE(16,*)
       WRITE(16,*)
       L=126.32D0/2
       CALL TRAPEZOIDS(-L,L,18,DENSITAT,RESULTATT)
       CALL SIMPSON(-L,L,18,DENSITAT,RESULTATS)
       WRITE(16,*) "Massa calculada per trapezis i Simpson:"
       WRITE(16,*) RESULTATT,RESULTATS
       CLOSE(16)

C 3)
C COMPARACIÓ DELS RESULTATS
C   a)
       OPEN(17,FILE="P4_1920_res2.dat")
       WRITE(17,*) "h, resultat trapezis i resultat Simpson (àrea)"
       DO K=4,22,1
        H=2.D0*PI/(2**K)
        CALL TRAPEZOIDS(-PI,PI,K,LONGITUD,RESULTATT)
        CALL SIMPSON(-PI,PI,K,LONGITUD,RESULTATS)
        WRITE(17,18) H,RESULTATT,RESULTATS
       ENDDO
18     FORMAT(3(E14.8,2X))
       CLOSE(17)

C   b)
       OPEN(19,FILE="P4_1920_res3.dat")
       WRITE(19,*) "h, resultat trapezis i resultat Simpson (massa)"
       DO K=4,22,1
        H=2.D0*L/(2**K)
        CALL TRAPEZOIDS(-L,L,K,DENSITAT,RESULTATT)
        CALL SIMPSON(-L,L,K,DENSITAT,RESULTATS)
        WRITE(19,18) H,RESULTATT,RESULTATS
       ENDDO
       CLOSE(19)


C ESTUDI DE L'ERROR
C   a)
       OPEN(20,FILE="P4_1920_res.dat")
       WRITE(20,*) "h, error trapezis i error Simpson (àrea)"
C l'aproximació més exacta l'agafem amb n=23
       CALL TRAPEZOIDS(-PI,PI,23,LONGITUD,EXACTEA)
       DO K=4,22,1
        H=2.D0*PI/(2**K)
        CALL TRAPEZOIDS(-PI,PI,K,LONGITUD,RESULTATT)
        CALL SIMPSON(-PI,PI,K,LONGITUD,RESULTATS)
C es guarda l'error: diferència entre cada valor i el que prenem com exacte
        WRITE(20,18) H,ABS(RESULTATT-EXACTEA),ABS(RESULTATS-EXACTEA)
       ENDDO

       WRITE(20,*)
       WRITE(20,*)

C   b)
       WRITE(20,*) "h, error trapezis i error Simpson (massa)"
C l'aproximació més exacta l'agafem amb n=23
       CALL TRAPEZOIDS(-L,L,23,DENSITAT,EXACTEM)
       DO K=4,22,1
        H=2.D0*L/(2**K)
        CALL TRAPEZOIDS(-L,L,K,DENSITAT,RESULTATT)
        CALL SIMPSON(-L,L,K,DENSITAT,RESULTATS)
C es guarda l'error: diferència entre cada valor i el que prenem com exacte
        WRITE(20,18) H,ABS(RESULTATT-EXACTEM),ABS(RESULTATS-EXACTEM)
       ENDDO
       CLOSE(20)

C creem les figures amb un mateix script
       CALL SYSTEM("gnuplot fig.gnu")



C 4)
       OPEN(21,FILE="P4_1920_res4.dat")

C guardem els resultats
       WRITE(21,*) "h, resultat trapezis i resultat Simpson (canvi var)"
C en aquest cas, t va de -PI/2 a PI/2
       DO K=4,20,1
        H=PI/(2**K)
        CALL TRAPEZOIDS(-PI/2,PI/2,K,F3,RESULTATT)
        CALL SIMPSON(-PI/2,PI/2,K,F3,RESULTATS)
        WRITE(21,18) H,RESULTATT,RESULTATS
       ENDDO

       WRITE(21,*)
       WRITE(21,*)

C guardem els errors
       WRITE(21,*) "h, error trapezis i error Simpson (canvi var)"
C l'aproximació més exacta l'agafem amb n=21
       CALL TRAPEZOIDS(-PI/2,PI/2,21,F3,EXACTEM)
       DO K=4,20,1
        H=PI/(2**K)
        CALL TRAPEZOIDS(-PI/2,PI/2,K,F3,RESULTATT)
        CALL SIMPSON(-PI/2,PI/2,K,F3,RESULTATS)
C es guarda l'error: diferència entre cada valor i el que prenem com exacte
        WRITE(21,18) H,ABS(RESULTATT-EXACTEM),ABS(RESULTATS-EXACTEM)
       ENDDO
       CLOSE(21)

       PRINT*, "Error 5 ordres de magnitud més petit amb el canvi"


       END PROGRAM

C SUBRUTINES:
C Regla dels trapezis:
       SUBROUTINE TRAPEZOIDS(X1,X2,K,FUNCI,INTEGRAL)
       IMPLICIT NONE
       INTEGER K,I
       DOUBLE PRECISION X1,X2,INTEGRAL,H,XK,FK
C separació entre particions constant, amplada intervals:
       H=(X2-X1)/(2**K)
C calculem el valor de la fucnió als extrems i sumem la seva contribució a la integral
       CALL FUNCI(X1,FK)
       INTEGRAL=0.5D0*FK
       CALL FUNCI(X2,FK)
       INTEGRAL=INTEGRAL+0.5D0*FK
C sumem les contribucions del punts intermitjos
       DO I=1,2**K-1,1
        XK=X1+I*H
        CALL FUNCI(XK,FK)
        INTEGRAL=INTEGRAL+FK
       ENDDO
C finalment només queda multiplicar per H per tenir la integral per trapezis. No ho faig fins al final perquè H és petita i així no perdem precisió treballant amb nombres petits.
       INTEGRAL=INTEGRAL*H
       RETURN
       END

C Regla de Simpson:
       SUBROUTINE SIMPSON(X1,X2,K,FUNCI,INTEGRAL)
       IMPLICIT NONE
       INTEGER K,I
       DOUBLE PRECISION X1,X2,INTEGRAL,H,XK,FK
C amplada dels intervals:
       H=(X2-X1)/(2**K)
C calculem el valor de la fucnió als extrems i sumem la seva contribució a la integral
       CALL FUNCI(X1,FK)
       INTEGRAL=FK
       CALL FUNCI(X2,FK)
       INTEGRAL=INTEGRAL+FK
C sumem les contribucions del punts intermitjos
       DO I=1,2**K-1,1
        XK=X1+I*H
        CALL FUNCI(XK,FK)
C per punts d'I parella es multiplica per 2
        IF (MOD(I,2).EQ.0) THEN
         INTEGRAL=INTEGRAL+2.D0*FK
C per I imparella es multiplica per 4
        ELSE
         INTEGRAL=INTEGRAL+4.D0*FK
        ENDIF
       ENDDO
C falta multiplicar les sumes anteriors per H/3:
       INTEGRAL=INTEGRAL*H/3
       RETURN
       END

C Funció de la longitud:
       SUBROUTINE LONGITUD(X,F)
       IMPLICIT NONE
       DOUBLE PRECISION X,F,PI,QUADRAT
       PI=4.D0*DATAN(1.D0)
       QUADRAT=((DCOS(X-2))**2*DEXP(-X**2+DSIN(X)))**2
       F=QUADRAT*DSQRT(PI-X)*0.35D0
       RETURN
       END

C Funció de la massa:
       SUBROUTINE DENSITAT(X,F)
       IMPLICIT NONE
       DOUBLE PRECISION X,F,L
       L=126.32D0/2
       F=0.42D-3*DSQRT(1-(X/L)**2)*(1-(X/L))*((X/L)**2+(X/L)+1)
       RETURN
       END

C Funció amb el canvi de variable x=Lsin(t)
       SUBROUTINE F3(T,F)
       IMPLICIT NONE
       DOUBLE PRECISION T,F
       F=0.42D-3*(DCOS(T))**2*(1-DSIN(T))*((DSIN(T))**2+DSIN(T)+1)
       RETURN
       END
