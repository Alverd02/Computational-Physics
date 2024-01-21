PROGRAM G2022

IMPLICIT NONE

DOUBLE PRECISION :: sigma,pi,a,b,cotasup,p,boxsize,f,resultat_f,error_f,int2,temp,eps
DOUBLE PRECISION, ALLOCATABLE :: distribucio_p(:),xhis(:),vhis(:),errhis(:)
INTEGER :: nbox,i,ierr,ndades,ISEED
COMMON/CONSTANTS/sigma,pi
EXTERNAL p,f,f2

ISEED=20470586
CALL SRAND(ISEED)

! 1
! a)

sigma = dsqrt(2.d0)
pi = 4.d0*atan(1.d0)

ndades = 100000
allocate(distribucio_p(ndades))
a = -4.d0*sigma
b = 4.d0*sigma
cotasup = 1.d0/(sigma*dsqrt(2.d0*pi))

CALL accepta(ndades,distribucio_p,a,b,cotasup,p)

nbox = 100
allocate(xhis(nbox),vhis(nbox),errhis(nbox))
CALL histograma(ndades, distribucio_p, a, b, nbox, xhis, vhis, errhis, boxsize, ierr)

OPEN(11,file="apartat1a.dat")

WRITE(11,*) "Histograma"

DO i=1,nbox

WRITE(11,*) xhis(i),vhis(i),errhis(i)

END DO

WRITE(11,"(/)")

CALL montecarlo(distribucio_p,ndades,f,P,resultat_f,error_f)

WRITE(11,*) "# Integral de f(x)"
WRITE(11,*) resultat_f,error_f
WRITE(11,"(/)")
WRITE(11,*) "# Convergencia"
    temp=0D0
    int2=10D0
    i=1
    eps = 1D-10

    !Mediante un do while vamos aumentando las iteraciones hasta obtener la precision deseada.
    !Para ello guardamos en la variable temporal "temp" el resultado de la anterior integral.
    !Si la diferencia entre las dos integrales es pequeña, el bucle acaba.

    do while (abs(temp-int2).GT. eps)
        temp = int2
        call simpsontresvuit(0D0,dsqrt(2D0),10*i,f2,int2)
        i = i+1 
        write(11,*) 100*i ,int2
    enddo
    write(11,*)
    write(11,*)
CLOSE(11)

END PROGRAM G2022

DOUBLE PRECISION FUNCTION f(x)

IMPLICIT NONE

DOUBLE PRECISION :: f,x

f = exp(-x**2)*(x**2-2.*x*dsin(x))

END FUNCTION

DOUBLE PRECISION FUNCTION f2(x)

IMPLICIT NONE

DOUBLE PRECISION :: f2,x

f2 = (x**2*dsin(x))/(x**3+dcos(x**2))

END FUNCTION

DOUBLE PRECISION FUNCTION p(x)

IMPLICIT NONE

DOUBLE PRECISION :: x,p,sigma,pi
COMMON/CONSTANTS/sigma,pi

p = exp(-x**2/(2.d0*sigma**2))/(sigma*dsqrt(2.d0*pi))

END FUNCTION

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

subroutine montecarlo(valors,n,func,dist,resultat,error)

IMPLICIT NONE

INTEGER :: n,i
DOUBLE PRECISION :: h,func,resultat,error,x2,dist
DOUBLE PRECISION, DIMENSION(n) :: valors

resultat = 0.d0
x2 = 0.d0

DO i=1,n 

h = func(valors(i))/dist(valors(i))
resultat = resultat + h
x2 = x2 + h**2

END DO

error = (1/dsqrt(dble(n)))*dsqrt((x2/dble(n)) - (resultat/dble(n))**2)
resultat = resultat/dble(n)

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