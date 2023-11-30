! © Marc I el Modest, 2019

!Important! Aquesta subrutina requereix la subrutina eqdif, que depèn del
!	problema, i subministra les derivades en funció de les altres coordenades.

subroutine Euler(ndim, yin, yout, x, h)
	implicit none
	integer ndim
	real*8 x, h
	real*8, dimension(ndim) :: yin, yout, dyin
!ndim: nombre de variables dependents del sistema.
!yin: conjunt de variables dependents.
!yout: conjunt de variables dependents després del pas.
!dyin: derivades de les variables dependents respecte la variable independent.
!x: variable independent del sistema.
!h: pas a la variable independent.
! Cridem la subrutina que ens dóna els valors per les variables dependents.
	call eqdif(ndim, yin, dyin, x)
! Realitzem el pas.
	yout = yin + h*dyin
end subroutine Euler

