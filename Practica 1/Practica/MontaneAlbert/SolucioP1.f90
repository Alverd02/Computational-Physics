!--------------------------------------------------------
!--              PRACTICA 1 
!--------------------------------------------------------
! Versió C1 curs 2018-19 
!
! ENTRADA:
! --------
! El programa demana el numero de linies 'k' de tall.                                                                                                                                   
!
! SORTIDA:
! --------
! Per pantalla dona informacions demanades en el guio
! de la practica.
! Genera un arxiu per tal de poder dibuixar grafiques.
!
! Programacio
! -----------
! 
! Llenguatge Fortran-90
! 
! Fa servir functions per calcular el numero 
! maxim de talls que es poden fer amb k linies rectes i 
! la suma demanada a l'enciat de la pràctica.
! 
!---------------------------------------------------------
! Barcelona 26-Sep-2018 R. Mayol
!---------------------------------------------------------

program practica1
implicit none 

!..............................
!... Definicio de variables
!..............................

integer (kind=4)  :: k,n,m    ! Indexs

!..............................
! Banner bonic del programa
!..............................

write(*,*) '           ****************'
write(*,*) '           ** PRACTICA 1 **'
write(*,*) '           **  Vers. C1  **'
write(*,*) '           **  2018-19   **'
write(*,*) '           ****************'

!-------------------------------------------------------------
!                      Apartat 1. 
!  -Es llegeix un valor de k entre 5 i 301 
!  -Calcula i escriu a la pantalla Pk
!-------------------------------------------------------------

k = 0
do while (k<5 .or. k > 301)
    write(*,*) 'Dona''m un valor entre 5 i 301'
    read(*,*,err=999,end=998) k 
end do 
998 continue

write(*,*) '   Amb k=  ',k,' el valor de Pk és= ',Pk(k)

!-------------------------------------------------------------
!                      Apartat 2. 
!  -Calcula pels valors M=14 i N=70 el valor del Sumatori  
!-------------------------------------------------------------

n = 70
m = 14 

write(*,*) '   Amb M= ',m,' i N= ',n,' el valor de S és= ',snm(n,m)
!----------------------------------------------------------------
!                      Apartat 3. 
!  -Construccio de l'arxiu P1-1819-res1.dat
!   dos columnes N, Sumatori des de M=3 fins N amb N=4,6,...,210
!----------------------------------------------------------------

open(7,file='P1-1819-res1.dat')
write(*,*) '   Generació de l''arxiu P1-1819-res1.dat'
write(7,*)'# Dues columnes: N, Suma(n,m) m=3'

m = 3 
do n=4,210,2
   write(7,*) N, snm(n,m)
end do

close(7)
write(*,*) '   Fi de programa '

!....................................
! S O R T I D A 
!....................................

stop 

999 stop 'Error a l''entrada de dades. El programa s''atura'

contains 

!...................................................
!   Funcio que calcula el valor de l'expressio Pk
!...................................................
! Valors de la constant Pi tret del llibre Schaum 
!
double precision function pk(k)
integer (kind=4) :: k
real    (kind=8) :: pi =3.1415926435897932384626d0
pk=(5.d0*k**2/3.d0)+pi-2.0d0*k
return 
end function pk

!.....................................
!...      Function Snm
!.....................................
!
double precision function Snm(n,m)
implicit none
real    (kind=8) :: p
integer (kind=4) :: n,m,k

Snm=0.0d0

do k=m,n
    snm = snm+Pk(k)
end do

return 
end function Snm

end program practica1


