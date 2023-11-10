PROGRAM P5

IMPLICIT NONE

INTEGER :: ISEED,ndades,nbox,i,ndat
DOUBLE PRECISION :: pi,e,a,b,func,cotasup, boxsize, ierr,xlam,valor_mitjag,valor_mitjat,varianciag,varianciat,desv_estg,desv_estt
DOUBLE PRECISION, dimension(50) :: xhis, vhis, errhis
DOUBLE PRECISION, dimension(110) :: xhis2, vhis2, errhis2
DOUBLE PRECISION, dimension(40000) :: numeros,xi
DOUBLE PRECISION, dimension(24000) :: xnumexpo,xi2
COMMON/CONSTANTS/pi,e
EXTERNAL func

ISEED=20470586
CALL SRAND(ISEED)

pi = 3.1415926535898 
e = 2.718281828459 

ndades = 40000
a = -pi
b = pi
cotasup = 0.4d0

OPEN(11,file="P5-23-24-res.dat")

CALL accepta(ndades,numeros,a,b,cotasup,func)
xi = numeros
nbox = 50
CALL histograma(ndades, xi, a, b, nbox, xhis, vhis, errhis, boxsize, ierr)
WRITE(11,"(/)")
WRITE(11,*)"# 1"
DO i=1,nbox
WRITE(11,*) xhis(i),vhis(i),errhis(i)
END DO

ndat = 24000
xlam = 6.d0/7
CALL sexponencial(ndat,xlam,xnumexpo)

valor_mitjag = sum(xnumexpo)/ndat
varianciag = sum((xnumexpo - valor_mitjag)**2)/ndat
desv_estg = sqrt(varianciag)

valor_mitjat = 1/xlam
varianciat = 1/xlam**2
desv_estt = sqrt(varianciat)

WRITE(11,"(/)")
WRITE(11,*) "# 2"
WRITE(11,*) valor_mitjag,valor_mitjat,varianciag,varianciat,desv_estg,desv_estt

a = 0
b = 6.d0*xlam
xi2 = xnumexpo
nbox = 110
CALL histograma(ndat, xi, a, b, nbox, xhis2, vhis2, errhis2, boxsize, ierr)
WRITE(11,"(/)")
WRITE(11,*)"# 3"
DO i=1,nbox
WRITE(11,*) xhis2(i),vhis2(i),errhis2(i)
END DO

CLOSE(11)

END PROGRAM P5

SUBROUTINE accepta(ndades,numeros,xlow,xhigh,cotasup,funcio)

IMPLICIT NONE

DOUBLE PRECISION :: xlow,xhigh,cotasup,funcio,p,valor_mitja,variancia,desv_est
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

valor_mitja = sum(numeros)/ndades
variancia = sum((numeros - valor_mitja)**2)/ndades
desv_est = sqrt(variancia)

WRITE(11,*)"# 0"
WRITE(11,*) valor_mitja,variancia,desv_est

RETURN
END


DOUBLE PRECISION FUNCTION func(x)

IMPLICIT NONE

COMMON/CONSTANTS/pi,e
DOUBLE PRECISION :: pi,e,x,y

y = (125*e**pi*x**2*(dsin(x))**2*e**(-abs(x)))/(4*(68*e**pi-70*pi-25*pi**2-68))

RETURN
END

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

SUBROUTINE sexponencial(ndades,xlam,xnumexpo)

IMPLICIT NONE

DOUBLE PRECISION,dimension(ndades) :: xnumexpo
INTEGER :: ndades,i
DOUBLE PRECISION :: xlam

DO i=1,ndades
xnumexpo(i) = (-1.d0*log(RAND()))/xlam
END DO

END