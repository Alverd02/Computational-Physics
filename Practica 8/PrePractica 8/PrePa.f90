program prepra8
    implicit none
    double precision :: E, V, pi, E1, E2, err, E3, Int
    integer :: i, N
    common/energies/E,V
    external derivades
    open(15,file = 'resultatsprepra8.dat')

    pi = 4.d0 * atan(1.d0)
    V = -2.4

    N = 400
    err = 10.d-5
    do i =  1,4
        E1 = (i*pi)**2/2.d0 + V +0.2d0
        E2 = (i*pi)**2/2.d0 + V +0.3d0
        call tir(E1, E2, N, err, E3)
        call normalitzar(E3, N)
    end do

    N = 20
    do i =  1,4
        E1 = (i*pi)**2/2.d0 + V +0.2d0
        E2 = (i*pi)**2/2.d0 + V +0.3d0
        call tir(E1, E2, N, err, E3)
        call normalitzar(E3, N)
    end do

    close(15)
end program prepra8



subroutine RungeKutta4(t0, dt, nequ, yyin, yyout,derivades)
    implicit none
    double precision :: t0, dt, t
    double precision, dimension(nequ) :: yyin, yyout, y, k1, k2, k3, k4
    integer :: nequ, i


    call derivades(nequ, t0, yyin, k1) !

    y = yyin + 0.5d0 * k1 * dt 

    t = t0 + 0.5d0*dt 
    call derivades(nequ, t0, y, k2)


   
        y = yyin + 0.5d0 * k2 * dt 
    t = t0 + 0.5d0*dt 
    call derivades(nequ, t0, y, k3)


  

        y = yyin + k3 * dt 
    t = t0 + dt 
    call derivades(nequ, t0, y, k4)


    !y(y0+h) = y0 + h/6*(k1 + 2*k2 + 2*k3 + k4)
        yyout = yyin + dt/6.d0*(k1+2.d0*k2+2.d0*k3+k4)


end subroutine

!Subroutina que calcula les derivades del vector yin i les retorna a dyout
!Pel cas particular amb l'equació d'Schrödinguer on y = (phi, dphi/dx)
!Amb notació dels apunts per RK4: t = x, yin = y, dyout = f(t,yin) = f(x,y)
subroutine derivades(nequ, t, yin, dyout)
    implicit none
    double precision :: t, yin(nequ), dyout(nequ), E, V
    integer :: nequ
    common/energies/E,V

    !Estem treballant en un cas on yin = [phi, dphi/dx] (mirar exemple pèndol apunts)
    !i volem dyout = [dphi/dx, d^2phi/dx^2]
    !La segona derivada ve donada per l'equació d'Schrödinguer

    dyout(1) = yin(2)
    dyout(2) = 2.d0 * (V-E) * yin(1)

end subroutine


!Subroutina que, donats uns límits d'integració, un nombre d'intervals
!i una energia E0, ens retorna una llista amb els phi, la phi final,
!I la Integral numèrica en aquest interval de phi^2dx
subroutine integrar(a, b, N, E0, Integral, phi_list, phi_f)
    implicit none
    double precision :: a, b, E0,  dx, x, Integral, E, V, phi_f
    double precision, dimension (2) :: funcin, phiRK
    double precision, dimension (N) :: phi2_list, phi_list
    integer :: N, Nequ, i, k
    common/energies/E,V
    external derivades


    funcin = [0.0d0,0.15d0]
    Nequ = 2
    E = E0
    dx = (b-a)/dble(N) !h

    do i = 1,N
        x = a + dx*i
        phi_list(i) = funcin(1)
        call RungeKutta4(x, dx, Nequ, funcin, phiRK,derivades) !Següent pas
        funcin = phiRK
    end do

    phi_f = phiRK(1)


end subroutine


subroutine tir(E1, E2, N, error, E3)
    implicit none
    double precision :: E1, E2, E3, error, xi, xf, phi1, phi2, phi3
    double precision :: I1, I2, I3, phi_list(N)
    integer :: N, iteracions, i

    xi = 0.d0
    xf = 1.d0
    iteracions = 20

    write(*,*) "Tir!"

    do i = 1, iteracions

        !Pas 2: Integrem
        call integrar(xi, xf, N, E1, I1, phi_list, phi1)
        call integrar(xi, xf, N, E2, I2, phi_list, phi2)

        !Pas 3: E3
        E3 = (E1*phi2-E2*phi1)/(phi2-phi1)
        
        !Pas 4: Estudiem si ha convergit i, si escau, repetim
        call integrar(xi, xf, N, E3, I3, phi_list, phi3)
        write(*,*) phi3
        if (abs(phi3) .lt. error) then
            write(*,*) "Ha convergit", i, E3
            exit
        else
            E1 = E2
            E2 = E3
        end if

    end do

end subroutine

!Subroutina que, a partir d'un interval (x1,x2), una llista i un enter N de forma
!que: iteracions = N , et retrona una aproximació numèrica del valor de la 
!intergral mitjançant el mètode de simpson amb 3 punts.
SUBROUTINE  simpson_list(x1,x2,k,func,resultat)

IMPLICIT NONE

DOUBLE PRECISION :: x1,x2,resultat,h,valor,x,valor0,valorN
double precision,dimension(k) :: func
INTEGER :: intervals,k,i

intervals = k
h = (x2-x1)/intervals

VALOR0 = func(1)
VALORn = func(k)

resultat = valor0 + valorN

DO i=1,intervals-1
x = x1
IF (MOD(i,3).eq.0) THEN

valor = func(i+1)
resultat=resultat+2*valor

ELSE

valor =  func(i+1)
resultat=resultat+3*valor

END IF

END DO

resultat = (3*h*resultat)/8.

RETURN
END 

subroutine trapezoids(x1,x2,N,y,integral)
    implicit none
    double precision x1,x2,x,integral,y(N),h
    integer N,int,i
    int = N
    integral = 0
    h = (x2-x1)/int
    do i=0,int
        x = x1+i*h
        integral = integral + h*y(i+1)
    enddo
    integral = integral - (y(1)+y(N))*h/2

end



subroutine normalitzar(E1, N)
    implicit none
    double precision :: phi2(N), llista_n(N),  Int, a, b, E1,phi(n)
    double precision :: x, dx, phi_F
    integer :: N, i

    a = 0.d0
    b = 1.d0

CALL integrar(a, b, N, E1, Int, phi, phi_f)

    dx = (b-a)/(1D0*N)
    do i = 1,N
        phi2(i) = phi(I)**2
    end do
    
    !Integrem
    call simpson_list(a, b, N, phi2, Int)
    !call trapezoids(a, b, N, phi2_list, Integral)

    write(15,*) "#E =", E1, "N =", N
    do i = 1,N
        x = a + dx*i
        llista_n(i) = PHI(i)/sqrt(Int)
        write(15,*) x,phi(i),llista_n(i)
    end do

write(15,"(/)")
end subroutine