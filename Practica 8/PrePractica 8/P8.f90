

program practica8


    implicit none
    double precision :: delta, coef, E1, E2, E3, E4, L, E, beta, xi, h
    double precision :: V, phi0, dphi0, x0, x1, phif, x, dx, V0, Int
    double precision :: a, b, E5, E6, error, Eout1, Eout2, Eout3, Eout4
    double precision :: beta1, beta2, beta3, Eout5, Eout6, Eout7, I4, I5, I6
    double precision, allocatable :: phi1(:), phi2(:), phi3(:), phi4(:)
    double precision, allocatable :: phi5(:), phi6(:), phi7(:)
    double precision, allocatable :: phin1(:), phin2(:), phin3(:)
    double precision, allocatable :: phin4(:), phin5(:), phin6(:)
    double precision, allocatable :: phin4quad(:), phin5quad(:), phin6quad(:)
    integer :: N, i, Neq
    common/dadesV/ delta, V0, xi
    common/energies/E,phi0,dphi0, coef
    external V

    open(15,file = "P8-22-23-res.dat")

    delta = 0.4d0 !A
    V0 = -50.d0 !eV
    coef = 3.80995d0 !eV A^2  (hbarra^2/m_e)
    L = 14.d0 !A

    phi0 = 0.d0 !A^(-1/2)
    dphi0 = 2.d-6 !A^(-3/2)

    x0 = -L/2.d0
    x1 = L/2.d0
    
    beta = 0.d0
    xi = 1.d0 ! A


    !______________Apartat 1______________________
    write(15,*) "#Apartat 1"
    N = 400
    allocate(phi1(N), phi2(N), phi3(N), phi4(N))

    E1 = -31.d0 !eV
    E2 = -30.d0 !eV
    E3 = -14.d0 !eV
    E4 = -13.d0 !eV

    Neq = 2

    call integrar(x0, x1, N, E1, Int, phi1, phif , V, beta)
    call integrar(x0, x1, N, E2, Int, phi2, phif , V, beta)
    call integrar(x0, x1, N, E3, Int, phi3, phif , V, beta)
    call integrar(x0, x1, N, E4, Int, phi4, phif , V, beta)

    dx = (x1-x0)/dble(N)
    do i=1,N
        x = x0 + i*dx
        write(15, *) x, phi1(i), phi2(i), phi3(i), phi4(i)
    end do

    write(15,*) ""
    write(15,*) ""
    write(15,*) ""


    !_________________Apartat 2______________________________
    !a)
    write(15,*) "#Apartat 2"

    E5 = -4.d0 !eV
    E6 = -3.5d0 !eV

    error = 1.d-6 !A^(-1/2)

    call tir(E1, E2, N, error, Eout1, x0, x1, V, beta)
    call tir(E3, E4, N, error, Eout2, x0, x1, V, beta)
    call tir(E5, E6, N, error, Eout3, x0, x1, V, beta)

    !b)
    allocate(phin1(N), phin2(N), phin3(N))

    call normalitza(x0, x1, Eout1, N, Int, phin1, V, beta)
    call normalitza(x0, x1, Eout2, N, Int, phin2, V, beta)
    call normalitza(x0, x1, Eout3, N, Int, phin3, V, beta)


    !_________________Apartat 3______________________________
    !a)

    allocate(phin4(N), phin5(N), phin6(N))
    beta1 = 0.d0 !eV
    beta2 = 5.d0 !eV
    beta3 = 15.d0 !eV
    write(*,*) ""
    write(*,*) ""
    write(*,*) "Apartat 3"

    write(*,*) "L'energia per la beta1 sera:"
    call tir(E1, E2, N, error, Eout5, x0, x1, V, beta1) 

    write(*,*) "L'energia per la beta2 sera:"
    call tir(E1, E2, N, error, Eout6, x0, x1, V, beta2)

    write(*,*) "L'energia per la beta3 sera:"
    call tir(E1, E2, N, error, Eout7, x0, x1, V, beta3)

    call normalitza(x0, x1, Eout5, N, Int, phin4, V, beta1)
    call normalitza(x0, x1, Eout6, N, Int, phin5, V, beta2)
    call normalitza(x0, x1, Eout7, N, Int, phin6, V, beta3)


    !b)

    !La probabilitat entre dos intervals a,b ve donada per la
    !Integral de phi^2 dx en a,b

    I4 = 0.d0
    I5 = 0.d0
    I6 = 0.d0
    h = (x1-x0)/N

    !Integrem fent servir trapezis
    do i = 1,N
        x = x0 + i*h
        if ((x .gt. -delta) .and. (x .lt. delta)) then
            I4 = I4 + phin4(i)**2*h
            I5 = I5 + phin5(i)**2*h
            I6 = I6 + phin6(i)**2*h
        end if
    end do

    open(16, file = "P8-22-23-res1.dat")
    write(16,*) "Beta = ", beta1, "Probabilitat =", I4
    write(16,*) "Beta = ", beta2, "Probabilitat =", I5
    write(16,*) "Beta = ", beta3, "Probabilitat =", I6

    close(16)
    close(15)
end program practica8




!___________________________________________________________________________________________
!________________________SUBROUTINES________________________________________________________



!Subtroutina que fa un pas de y(t0) a y(t0+dt) amb el mètode de RK4
!Amb notació dels apunts per a RK4: t0 = x0, dt = h, yyin = y0, yout = y(x0+h)
subroutine RungeKutta4(t0, dt, nequ, yyin, yyout, V, beta)
    implicit none
    double precision :: t0, dt, t, V, beta
    double precision, dimension(nequ) :: yyin, yyout, y, k1, k2, k3, k4
    integer :: nequ, i

    !k1:
    call derivades(nequ, t0, yyin, k1, V, beta) !k1 = f(x0,y0)


    !k2 = f(x0 + h/2, y0 + h/2*k2)
    do i = 1, nequ 
        y(i) = yyin(i) + 0.5d0 * k1(i) * dt !Actualitzem el valor de les y
    end do
    t = t0 + 0.5d0*dt !Actualitzem el valor de les x
    call derivades(nequ, t0, y, k2, V, beta)


    !k3 = f(x0 + h/2, y0 + h/2*k2)
    do i = 1, nequ 
        y(i) = yyin(i) + 0.5d0 * k2(i) * dt !Actualitzem el valor de les y
    end do
    t = t0 + 0.5d0*dt !Actualitzem el valor de les x
    call derivades(nequ, t0, y, k3, V, beta)


    !k4 = f(x0 + h, y0 + h*k3)
    do i = 1, nequ 
        y(i) = yyin(i) + k3(i) * dt !Actualitzem el valor de les y
    end do
    t = t0 + dt !Actualitzem el valor de les x
    call derivades(nequ, t0, y, k4, V, beta)


    !y(y0+h) = y0 + h/6*(k1 + 2*k2 + 2*k3 + k4)
    do i = 1, nequ
        yyout(i) = yyin(i) + dt/6.d0*(k1(i)+2.d0*k2(i)+2.d0*k3(i)+k4(i))
    end do

end subroutine




!Subroutina que calcula les derivades del vector yin i les retorna a dyout
!Pel cas particular amb l'equació d'Schrödinguer on y = (phi, dphi/dx)
!Amb notació dels apunts per RK4: t = x, yin = y, dyout = f(t,yin) = f(x,y)
subroutine derivades(nequ, t, yin, dyout, V, beta)
    implicit none
    double precision :: t, yin(nequ), dyout(nequ), E, V, phi0, dphi0, coef, beta
    integer :: nequ
    common/energies/E,phi0,dphi0,coef

    !Estem treballant en un cas on yin = [phi, dphi/dx] (mirar exemple pèndol apunts)
    !i volem dyout = [dphi/dx, d^2phi/dx^2]
    !La segona derivada ve donada per l'equació d'Schrödinguer

    dyout(1) = yin(2)
    dyout(2) =  (V(t, beta)-E) * yin(1)/coef

end subroutine





!Subroutina que, donats uns límits d'integració, un nombre d'intervals
!i una energia E0, ens retorna una llista amb els phi, la phi final,
!I la Integral numèrica en aquest interval de phi^2dx
subroutine integrar(a, b, N, E0, Integral, phi_list, phi_f, V, beta)
    implicit none
    double precision :: a, b, E0,  dx, x, Integral, E, V, phi_f
    double precision :: phi0, dphi0, coef, beta
    double precision, dimension (2) :: phi, phiRK
    double precision, dimension (N) :: phi2_list, phi_list
    integer :: N, Nequ, i, k

    common/energies/E,phi0,dphi0, coef

    !Valors inicials
    phi(1) = phi0
    phi(2) = dphi0

    Nequ = 2
    E = E0
    dx = (b-a)/dble(N) !h

    do i = 1,N
        x = a + dx*i
        phi2_list(i) = phi(1)**2
        phi_list(i) = phi(1)
        call RungeKutta4(x, dx, Nequ, phi, phiRK, V, beta) !Següent pas
        do k = 1, Nequ
            phi(k) = phiRK(k)
        end do
    end do

    phi_f = phiRK(1)
    
    !Integrem
    call simpson_list(a, b, N, phi2_list, Integral)

end subroutine




!Subroutina que troba els autovalors energia de l'equació d'Schrödinguer
!Amb el mètode de tir amb les condicions donades
subroutine tir(E1_, E2_, N, error, E3, a, b, V, beta)
    implicit none
    double precision :: E1_, E2_, E1, E2, E3, error, xi, xf, phi1, phi2, phi3
    double precision :: I1, I2, I3, phi_list(N), a, b, V, beta
    integer :: N, iteracions, i


    iteracions = 500

    E1 = E1_
    E2 = E2_
    
    write(*,*) "Tir!"

    do i = 1, iteracions

        !Pas 2: Integrem
        call integrar(a, b, N, E1, I1, phi_list, phi1, V, beta)
        call integrar(a, b, N, E2, I2, phi_list, phi2, V, beta)

        !Pas 3: E3
        E3 = (E1*phi2-E2*phi1)/(phi2-phi1)
        
        !Pas 4: Estudiem si ha convergit i, si escau, repetim
        call integrar(a, b, N, E3, I3, phi_list, phi3, V, beta)
        write(15,*) i,E3
        if (abs(phi3) .lt. error) then
            write(*,*) "Ha convergit", i, "Energia (eV) :", E3
            exit
        else
            E1 = E2
            E2 = E3
        end if

    end do

    write(15,*) ""
    write(15,*) ""
    write(15,*) ""


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




!Subroutina per a normalitzar l'equació d'schrödinguer a una energia
!i límits d'integració donats
subroutine normalitza(a, b, E1, N, Int, llista_n, V, beta)
    implicit none
    double precision :: llista(N), llista_n(N),  Int, a, b, E1
    double precision :: x, dx, phi, V, beta
    integer :: N, i
    
    dx = (b-a)/(1D0*N)

    call integrar(a, b, N, E1, Int, llista, phi, V, beta)
    write(15,*) "#E =", E1, "N =", N
    do i = 1,N
        x = a + dx*i
        llista_n(i) = llista(i)/sqrt(Int)
        write(15,*) x,llista(i),llista_n(i)
    end do
    write(15,*) ""
    write(15,*) ""
end subroutine






double precision function V(x, beta)
    implicit none
    double precision :: x, delta, V0, xi, beta
    common/dadesV/ delta, V0, xi
    V = V0 * sinh(2.d0) / (cosh(2.d0) + cosh(x/delta)) + beta*sin(x/xi)
end function