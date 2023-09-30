C PROGRAMA (MUY) SENCILLO PARA PROBAR EL COMPILADOR
C
C BJD SEP 2015
C
C LAST MODIFIED: 21 SEP 2015
C
      IMPLICIT NONE
      REAL X1,X2,X3

        X1=2.0E30
        X2=2.0E8
        X3=X1+X2
        WRITE(*,*) ' Redondeo'
        WRITE(*,*) ' NUMEROS    x1=',x1,' x2=',x2
        WRITE(*,*) ' X3-X1=X2=',X3-X1
        WRITE(*,*) ' X2=',X2
        WRITE(*,*)
        WRITE(*,*) ' PORQUE EN UN CASO X2=0 Y EN EL OTRO NO?'
        END
