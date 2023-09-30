C PROGRAM BISECCION
C
C THE PROGRAM LOOKS FOR A SOLUTION TO F(X)=0. 
C THE FUNCTION IS DEFINED EXTERNAL TO THE ROUTINE
C THE TWO INITIAL VALUES, A, B SHOULD FULFILL F(A) F(B)<0
C THE PROGRAM STOPS EITHER IF THE EXACT POINT IS FOUND OR
C IF THE DESIRED ACCURACY HAS BEEN ACHIEVED

C BJD SEP 2015
C LAST REVISED 4 OCT 2015
C
       IMPLICIT NONE
c extremos del intervalo y punto central
       REAL A,B,C,F

c precision requerida y error
       REAL EPS,DIFF
c un contador y un numero maximo de iteraciones
       INTEGER I,MAXITER


C PRECISION REQUERIDA
       EPS=0.0001

C VALORES EXTREMOS INICIALES
       A=3.
       B=8.
       WRITE(*,500) A,B
500    FORMAT('A=',F9.3,2X,'B=',F9.3)
      
C CALCULA MAXITER (le sumo uno para asegurar)

	MAXITER=NINT(LOG((B-A)/EPS)/LOG(2.))+1
        WRITE(*,*) "MAXITER=",MAXITER
C COMPRUEBA F(A)F(B)<0

          IF (F(A)*F(B).GE.0.) THEN 
            PRINT*,"LA FUNCION NO CAMBIA DE SIGNO EN A,B",
     1   A,B,F(A),F(B)
            STOP
          ENDIF

C COMIENZA EL METODO
       DO I=1,MAXITER
          C=(A+B)/2.        

          IF (F(C).eq.0.) THEN 
           PRINT*,"SOLUCION EXACTA X=",C
           STOP
          ENDIF
          IF (F(A)*F(C).LT.0.) THEN 
             B=C
           ELSE
             A=C          
          ENDIF
          DIFF=(B-A)

      IF (DIFF.LE.EPS) THEN 
         PRINT*,"SOLUCION APROXIMADA X=",C
         PRINT*,"ERROR <",DIFF
         STOP
      ENDIF

       PRINT*,"ITERACION NUMERO:",I,"C=",C, " error=",diff 
	ENDDO
       END

c funcion utilizada
       FUNCTION F(X)
       IMPLICIT NONE
       REAL X,F
       F=SIN(X)**2*X-1.
       END
