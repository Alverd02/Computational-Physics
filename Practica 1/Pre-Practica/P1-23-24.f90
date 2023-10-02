PROGRAM P1

IMPLICIT NONE

REAL :: E1 = 3.72
REAL :: E_fermi = 0.
INTEGER :: k,i,N
LOGICAL :: isNotValid = .TRUE.

! Apartat 1

DO WHILE (isNotValid)
    WRITE(*,*) "Escriu un numero entre 2 i 40"
    READ(*,*) k
    IF (k.GT.40.OR.k.LE.2) THEN
        WRITE(*,*) "Numero fora del rang. Escriu un altre."
    ELSE
        WRITE(*,*) "El valor Ek es: ", k**2*E1, "eV"
        isNotValid = .FALSE.
    ENDIF
END DO

! Apartat 2

DO i = 1,40
    E_fermi = E_fermi + i**2*E1
END DO
WRITE(*,*) "Energia de fermi: ", E_fermi, "eV"

! Apartat 3

OPEN(10,file = "P1-23-24.dat")
DO N=1,40
    E_fermi = 0
    DO i =1,N
        E_fermi = E_fermi + i**2*E1
    END DO
    WRITE(10,*) N,E_fermi
END DO

END PROGRAM P1