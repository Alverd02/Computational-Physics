PROGRAM P1

IMPLICIT NONE

INTEGER :: k,i,j
REAL :: S_asim, P_k = 0., n_0 = 28, n_1 = 65,e
COMMON/CONSTANTS/e 
e = 2.71828

! Apartat 1

WRITE(*,*) "Escriu un numero enter entre 15 i 221: "
READ(*,*) k
WRITE(*,*) "La funcio val: ", 3*k**2/5 + e + 10*k

! Apartat 2


DO i=28,65
    P_k = P_k + 3*i**2/5. + e + 10*i
END DO

WRITE(*,*) "El valor de S^M_N val (sense subrutina): ",P_k

CALL S(P_k, n_0, n_1)
WRITE(*,*) "El valor de S^M_N val (amb subrutina): ",P_k

! Apartat 3

OPEN(11,file = "P1-23-24-res1.dat")

DO i=8,311,3
    P_k = 0.
    DO j=8,i
        P_k = P_k + 3*j**2/5. + e + 10*j
    END DO
    WRITE(11,*) i,P_k
END DO

! Apartat 5

OPEN(12,file = "P1-23-24-res2.dat")

DO i=8,311,3
    P_k = 0.
    S_asim = i**3/5
    DO j=8,i
        P_k = P_k + 3*j**2/5 + e + 10*j
    END DO
    WRITE(12,*) i,P_k/S_asim
END DO

END PROGRAM P1

SUBROUTINE S(P_k,n_0,n_1)
    COMMON/CONSTANTS/e 
    INTEGER :: n_0,n_1,i
    REAL :: P_k
    P_k = 0.
    DO i=n_0,n_1
        P_k = P_k + 3*i**2/5 + e + 10*i
    END DO 
    RETURN
END

