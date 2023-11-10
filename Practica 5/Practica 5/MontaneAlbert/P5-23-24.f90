PROGRAM P5

IMPLICIT NONE

INTEGER :: ISEED,nbox,i,ndades,k
DOUBLE PRECISION :: pi,e,L,cotasup,a,b,boxsize,ierr,x2,resultat,resultat2,sigma,mu,valor_mitjag,varianciag,desv_estg
DOUBLE PRECISION, dimension(50000) :: numeros,xi
DOUBLE PRECISION, dimension(20000) :: yres,xi2
DOUBLE PRECISION, dimension(100) :: xhis, vhis, errhis
COMMON/CONSTANTS/pi,e,L
EXTERNAL func

ISEED=20470586
CALL SRAND(ISEED)

pi = 3.1415926535898 
e = 2.718281828459 

! Apartat 1
! a)

L = 4.d0

ndades = 50000
a = -L*pi
b = L*pi
cotasup = 1/(L*pi)

OPEN(11,file="P5-23-24-res.dat")

CALL accepta(ndades,numeros,a,b,cotasup,func)
xi = numeros
nbox = 100
CALL histograma(ndades, xi, a, b, nbox, xhis, vhis, errhis, boxsize, ierr)

WRITE(11,*) "# punts, valors i errors"
DO i=1,nbox
WRITE(11,*) xhis(i),vhis(i),errhis(i)
END DO
WRITE(11,"(/)")

WRITE(11,*) "# Probabilitat i comprobació normalització"

! b)

x2 = (L*pi)/2
k = 12

CALL simpsontresvuit(a,x2,k,func,resultat)
CALL simpsontresvuit(a,b,k,func,resultat2)
WRITE(11,*) resultat, resultat2
WRITE(11,"(/)")

! Apartat 2
! a)
ndades = 20000
sigma = 3.d0
mu = 0.d0
CALL boxmuller(ndades,sigma,mu,yres)
a = -4*sigma
b = 4*sigma
xi2 = yres
nbox = 100
CALL histograma(ndades, xi2, a, b, nbox, xhis, vhis, errhis, boxsize, ierr)

WRITE(11,*) "# punts, valors i errors"
DO i=1,nbox
WRITE(11,*) xhis(i),vhis(i),errhis(i)
END DO

! b)

valor_mitjag = sum(yres)/ndades
varianciag = sum((yres - valor_mitjag)**2)/ndades
desv_estg = sqrt(varianciag)


WRITE(11,"(/)")
WRITE(11,*) "# Comparacio valor mig, variancia i desviació"
WRITE(11,*) valor_mitjag,mu,varianciag,sigma**2,desv_estg,sigma

CLOSE(11)


END PROGRAM P5

! Subrutina extreta del campus virtual

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

DOUBLE PRECISION FUNCTION func(x)

IMPLICIT NONE

COMMON/CONSTANTS/pi,e,L
DOUBLE PRECISION :: pi,x,y,L,e

y = (dsin(x/L)**2)/(L*pi)

RETURN
END

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