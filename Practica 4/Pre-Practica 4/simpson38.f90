! © Marc I el Modest, 2019

subroutine simpson38 (x1, x2, k, fun, integral)
	implicit none
	integer i, k, k_
	real*8 x1, x2, integral, fun
!i: comptador.
!k: nombre d'intervals.
!k_: nombre d'intervals conveninentment fet múltiple de 3.
!x1, x2: límits d'integració.
!integral: doncs això.
!fun: funció a integrar.
!Assegurem que el nombre sigui parell mitjançant aritmètica d'enters.
	k_ = ((k+2)/3)*3
!Sumem els punts exteriors.
	integral = fun(x1) + fun(x2)
!Sumem els punts del tipus 3n + 1 amb el seu pes.
	do i=1, (k_ - 2), 3
		integral = integral + 3.d0*fun(x1 + (x2-x1)*dble(i)/dble(k_))
	enddo
!Sumem els punts del tipus 3n + 2 amb el seu pes.
	do i=2, (k_ - 1), 3
		integral = integral + 3.d0*fun(x1 + (x2-x1)*dble(i)/dble(k_))
	enddo
!Sumem els punts múltiples de 3, amb el seu pes.
	do i=3, (k_ - 3), 3
		integral = integral + 2.d0*fun(x1 + (x2-x1)*dble(i)/dble(k_))
	enddo
!Multipliquem el resultat per la longitud dels intervals i els paràmetres de Simpson -3/8.
	integral = 3.d0*(x2-x1)*integral/(8.d0*dble(k_))
end subroutine simpson38
