
program prepra8
    implicit none
    double precision :: E, V, pi, E1, E2, err, E3, Int
    double precision, allocatable :: phi_list(:)
    integer :: i, N
    common/energies/E,V

    open(15,file = 'resultatsprepra8.dat')

    pi = 4.d0 * atan(1.d0)
    V = -2.4

    N = 400
    allocate(phi_list(N))
    err = 10.d-5
    do i =  1,4
        E1 = (i*pi)**2/2.d0 + V +0.2d0
        E2 = (i*pi)**2/2.d0 + V +0.3d0
        call tir(E1, E2, N, err, E3)
        call normalitza(E3, N, Int, phi_list)
    end do

    close(15)
end program prepra8


!Amb notació dels apunts per a RK4: t0 = x0, dt = h, yyin = y0, yout = y(x0+h)
subroutine RungeKutta4(t0, dt, nequ, yyin, yyout)
    implicit none
    double precision :: t0, dt, t
    double precision, dimension(nequ) :: yyin, yyout, y, k1, k2, k3, k4
    integer :: nequ, i

    !k1:
    call derivades(nequ, t0, yyin, k1) !k1 = f(x0,y0)


    !k2 = f(x0 + h/2, y0 + h/2*k2)
    do i = 1, nequ 
        y(i) = yyin(i) + 0.5d0 * k1(i) * dt !Actualitzem el valor de les y
    end do
    t = t0 + 0.5d0*dt !Actualitzem el valor de les x
    call derivades(nequ, t0, y, k2)


    !k3 = f(x0 + h/2, y0 + h/2*k2)
    do i = 1, nequ 
        y(i) = yyin(i) + 0.5d0 * k2(i) * dt !Actualitzem el valor de les y
    end do
    t = t0 + 0.5d0*dt !Actualitzem el valor de les x
    call derivades(nequ, t0, y, k3)


    !k4 = f(x0 + h, y0 + h*k3)
    do i = 1, nequ 
        y(i) = yyin(i) + k3(i) * dt !Actualitzem el valor de les y
    end do
    t = t0 + dt !Actualitzem el valor de les x
    call derivades(nequ, t0, y, k4)


    !y(y0+h) = y0 + h/6*(k1 + 2*k2 + 2*k3 + k4)
    do i = 1, nequ
        yyout(i) = yyin(i) + dt/6.d0*(k1(i)+2.d0*k2(i)+2.d0*k3(i)+k4(i))
    end do

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
    double precision, dimension (2) :: phi, phiRK
    double precision, dimension (N) :: phi2_list, phi_list
    integer :: N, Nequ, i, k

    common/energies/E,V

    !Valors inicials
    phi(1) = 0.d0
    phi(2) = 0.15d0

    Nequ = 2
    E = E0
    dx = (b-a)/dble(N) !h

    do i = 1,N
        x = a + dx*i
        phi2_list(i) = phi(1)**2
        phi_list(i) = phi(1)
        call RungeKutta4(x, dx, Nequ, phi, phiRK) !Següent pas
        do k = 1, Nequ
            phi(k) = phiRK(k)
        end do
    end do

    phi_f = phiRK(1)
    
    !Integrem
    call simpson_list(a, b, N, phi2_list, Integral)
    !call trapezoids(a, b, N, phi2_list, Integral)

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
subroutine simpson_list(x1, x2, N, y_list, integral)
	implicit none
	double precision :: x1, x2, integral, x_i, x_m, x_f, delta_x, funci, h
    double precision, dimension(N) :: y_list
	integer :: k, i, intervals, N
	intervals = N
	h = (x2-x1)/intervals
	integral = 0.d0
    !El mètode diferencia tres intervals:
	do i = 0, intervals-1
        if((i .eq. 0) .or. (i .eq. intervals)) then
            integral = integral + y_list(i+1)* h/3.d0
        else if (mod(i,2) .eq. 0) then
            integral = integral + 2.d0*y_list(i+1)* h/3.d0 
        else if (mod(i,2) .eq. 1) then
            integral = integral + 4.d0*y_list(i+1)* h/3.d0
        end if
     end do
end subroutine

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



subroutine normalitza(E1, N, Int, llista_n)
    implicit none
    double precision :: llista(N), llista_n(N),  Int, a, b, E1
    double precision :: x, dx, phi
    integer :: N, i

    a = 0.d0
    b = 1.d0
    
    dx = (b-a)/(1D0*N)

    call integrar(a, b, N, E1, Int, llista, phi)
    write(15,*) "#E =", E1, "N =", N
    do i = 1,N
        x = a + dx*i
        llista_n(i) = llista(i)/sqrt(Int)
        write(15,*) x,llista(i),llista_n(i)
    end do
    write(15,*) ""
    write(15,*) ""
end subroutine