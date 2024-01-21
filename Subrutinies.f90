! Trobar 0 de la funció fa falta el valor de la funció i la seva derivada

SUBROUTINE newtonraphson(fun,xini,preci,nitera,valorarrel)
IMPLICIT NONE
DOUBLE PRECISION :: xini,preci,valorarrel,x1,valordf,valorf,eps,t
INTEGER :: nitera,i


DO i=1,999999

CALL fun(xini,valorf,eps,valordf,t)   
x1 = xini - valorf/valordf
IF (ABS(valorf/valordf).LE.preci) THEN
valorarrel = x1
nitera = i

RETURN 
STOP
END IF
valorarrel = x1
nitera = i
xini = x1
END DO
END

! Metode de biseccio per trobar 0  de funcions, nomes va falta el valor de la funció

SUBROUTINE biseccio(fun,A,B,preci,nitera,valorarrel)
IMPLICIT NONE
DOUBLE PRECISION :: A,B,preci,valorarrel,valorf_A,valorf_B,d,valordf_A,valordf_B,valorf_d,valordf_d
INTEGER :: i,max_nitera,nitera

max_nitera=INT(LOG((B-A)/preci)/LOG(2.0d0))+1

IF (valorf_A*valorf_B.GT.0.0d0) THEN
WRITE(*,*) "No se puede aplicar el metodo, no hay cambio de singo"
STOP
END IF

DO i=1,max_nitera
d = (A+B)/2.0d0
CALL fun(A,valorf_A,valordf_A)
CALL fun(B,valorf_B,valordf_B)
CALL fun(d,valorf_d,valordf_d)

IF (valorf_d.EQ.0.0d0) THEN
valorarrel = d
nitera = i
RETURN 
STOP
END IF

IF ((valorf_A*valorf_d).LT.0.0D0) THEN
B = d
ELSE 
A = d
END IF


IF ((B-A).LT.preci) THEN
valorarrel = d
nitera = i
RETURN 
STOP
END IF

END DO

END


! Meteode de trapezis

SUBROUTINE trapezoidalrule(x1,x2,k,func,resultat)

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,x,valor,resultat2,resultat1,valor_0,valor_N,func
INTEGER :: k,intervals,i

intervals = 3.d0**k
h = (x2-x1)/intervals

x=x1
valor_0 = func(x)
x = x2
valor_N =  func(x)

resultat1 = ((valor_0+valor_N)*h)/2.
resultat2 = 0.d0

DO i=1,(intervals-1)
x = x1
x = x+i*h
valor =  func(x)
resultat2 = valor + resultat2
END DO

resultat = resultat1 + h*resultat2
RETURN

END

! Calcula el valor d'una integral, d'una funció (retorna un valor), tambe es por arreglar perque la funció d'entrada sigui una llita de valors.


SUBROUTINE  simpsontresvuit(x1,x2,k,func,resultat)

! Metode de Simpsin 3/8 repetit, utilitzem 3^k intervals i começant des de 1 fins intervals-1 
!multipliquem per 2  o per 3 segons si i es dsivisible per 3 o no.

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,valor,x,valor0,valorN,func
INTEGER :: intervals,k,i

intervals = 3.d0**k
h = (x2-x1)/intervals

valor0 = func(x1)
valorN = func(x2)

resultat = valor0 + valorN

DO i=1,intervals-1
x = x1
IF (MOD(i,3).eq.0) THEN

valor =  func(x+i*h)
resultat=resultat+2*valor

ELSE

valor =  func(x+i*h)
resultat=resultat+3*valor

END IF

END DO

resultat = (3*h*resultat)/8.

RETURN
END 

!Genera uns histograma

subroutine histograma(ndat, xi, a, b, nbox, xhis, vhis, errhis, boxsize, ierr)
	implicit none
	integer :: ndat, nbox, ierr, i, ibox, icount
	double precision :: a, b, boxsize
	double precision, dimension(ndat) :: xi
	double precision, dimension(nbox) :: xhis, vhis, errhis

	if (a.ge.b) then 
		ierr=1
		return
	endif
	boxsize=(b-a)/nbox

	icount=0

	do i=1,nbox
		vhis(i)=0
		errhis(i)=0
	enddo

	do i=1,ndat
		if (xi(i).ge.a.and.xi(i).le.b) then 
			ibox=int((xi(i)-a)/boxsize)+1
		if (ibox.eq.nbox+1) ibox=nbox 

			vhis(ibox)=vhis(ibox)+1
			icount=icount+1
		endif
	enddo

	if (icount.eq.0) then 
		ierr=2
		return
	endif

	ierr=0
	print*,"accepted:",icount," out of:",ndat

	do i=1,nbox
		xhis(i)=a+boxsize/2.d0+(i-1)*boxsize
		errhis(i)=sqrt(vhis(i)/icount*(1.d0-vhis(i)/icount))/boxsize / sqrt(dble(icount))
		vhis(i)=vhis(i)/icount/boxsize
	enddo
end subroutine histograma

! Genera numeros aleatoris segons una distribució

SUBROUTINE accepta(ndades,numeros,xlow,xhigh,cotasup,funcio)

IMPLICIT NONE

DOUBLE PRECISION :: xlow,xhigh,cotasup,funcio,p
DOUBLE PRECISION, dimension(ndades) :: numeros
INTEGER :: ndades,n,i

n = 1

DO WHILE (n.lt.ndades)
numeros(n) = (xhigh-xlow)*RAND() + xlow
p = cotasup*RAND()
IF (funcio(numeros(n)).ge.p) then

n = n + 1

END IF
END DO

RETURN
END

! Genera numeros aleatoris gaussians

SUBROUTINE boxmuller(ndades,sigma,mu,yres)

IMPLICIT NONE
COMMON/CONSTANTS/pi
DOUBLE PRECISION,dimension(ndades) :: yres
INTEGER :: ndades,i
DOUBLE PRECISION :: mu,sigma,pi,x1,x2

DO i=1,ndades
x1 = RAND()
x2 = RAND()
yres(i) = mu + sigma*dsqrt(-2*dlog(x1))*dcos(2*pi*x2)
END DO

END

! Montecarlo cru vol dir que li pasem una funcio, els seus limits i retorna un valor i un error

SUBROUTINE montecarlocru(a,b,n,func,resultat,error)

IMPLICIT NONE

INTEGER :: n,i
DOUBLE PRECISION :: h,a,b,func,resultat,error,x2

resultat = 0.d0
x2 = 0.d0

DO i=1,n 

h = (b-a)*func((b-a)*rand()+a)
resultat = resultat + h
x2 = x2 + h**2

END DO

error = (1/dsqrt(dble(n)))*dsqrt((x2/dble(n)) - (resultat/dble(n))**2)
resultat = resultat/dble(n)

RETURN
END

! Integral de montecarlo: si la distribucio de "valors" es d'una desitibució dins de f(x) el meteode es queda com està, si f(x) es general i tenim valors d'una distribució externa
! h se l'haura de dividir per la distribució en el valor que toqui.( retorna un valor i el seru error)

subroutine montecarlo(valors,n,func,resultat,error)

IMPLICIT NONE

INTEGER :: n,i
DOUBLE PRECISION :: h,func,resultat,error,x2
DOUBLE PRECISION, DIMENSION(n) :: valors

resultat = 0.d0
x2 = 0.d0

DO i=1,n 

h = func(valors(i))
resultat = resultat + h
x2 = x2 + h**2

END DO

error = (1/dsqrt(dble(n)))*dsqrt((x2/dble(n)) - (resultat/dble(n))**2)
resultat = resultat/dble(n)

RETURN
END

! Runge-Kutta 4 resol un pas d'una equació diferencial d'un ordre de magnitud igual o mes que 2,

SUBROUTINE RungeKutta4(x0,dx,funcin,dfuncout,nequs,edofuncio)
IMPLICIT NONE
DOUBLE PRECISION :: dx,x,x0
DOUBLE PRECISION,DIMENSION(nequs) :: funcin, dfuncout, y, k1, k2, k3, k4
INTEGER :: nequs
CALL edofuncio(nequs, x0, funcin, k1) !
y = funcin + 0.5d0 * k1 * dx 
x = x0 + 0.5d0*dx 
CALL edofuncio(nequs, x, y, k2)
y = funcin + 0.5d0 * k2 * dx
x = x0 +  0.5d0*dx
CALL edofuncio(nequs, x, y, k3)
y = funcin + k3 * dx 
x = x0 + dx 
CALL edofuncio(nequs, x, y, k4)
dfuncout = funcin + dx/6.d0*(k1+2.d0*k2+2.d0*k3+k4)

END

SUBROUTINE EDO(nequs, x,funcin , dyoutput)
IMPLICIT NONE

DOUBLE PRECISION :: 
DOUBLE PRECISION ,DIMENSION(nequs)::funcin,dyoutput
INTEGER :: nequs
 

dyoutput(1) = funcin(2)
dyoutput(2) =  ! equació diferecial


END

SUBROUTINE integralRK4(x0,xf, N, E0, phi, phi_f)
IMPLICIT NONE
DOUBLE PRECISION :: x0,xf,E0,dx,x,E,V0,phi_f,delta,const
DOUBLE PRECISION,DIMENSION(2) :: funcin,phiRK
DOUBLE PRECISION,DIMENSION(N) ::  phi
INTEGER :: N, nequs, i
COMMON/CONSTANTS/E,V0,delta,const
EXTERNAL EDO
funcin = [0.0d0,2*1.0D-6]
nequs = 2
E = E0
dx = (xf-x0)/dble(N) 

DO i =1,n
x = x0 + dx*i
phi(i) = funcin(1)
call RungeKutta4(x,dx,funcin,phiRK,nequs,EDO)
funcin = phiRK
END DO
phi_f = phiRK(1)
END