program practica6

    implicit none
    double precision:: a, b, fun1, I1, s1, pi, L, M, xnums(1000000), p, fun2, i2, s2
    integer:: i, ndat
    external:: fun1, p, fun2
    common/dades/pi, L
    common/vectors/xnums

    !Apartat 1-a
    !Definim la variable 'pi'
    pi = 4*atan(1.d0)
    
    !Obrim el document on escriurem les dades
    open(1, file='P6-22-23-res.dat')
    write(1,'(A, /)') '#Apartat 1-a'

    !Definim els límits de l'integral
    a = 0
    b = 2*pi

    !Creem un bucle per escriure els resultats de cridar la subrutina 'montecarlo'
    do i = 10000, 1000000, 10000
        call montecarlo(a, b, fun1, i, I1, s1)
        write(1,*) i, I1, s1
    end do

    !Apartat 1-b
    !Definim les variables necessàries
    L = 50.d0
    M = 1/L
    a = 0.d0
    b = 2*L
    ndat = 1000000
    !Cridem la subrutina 'acceptrebuig' per calcular els 1000000 valors
    call acceptrebuig(ndat, xnums, a, b, M, p)

    !Apartat 1-c
    write(1,'(/,A,/)') '#Apartat 1-c'

    !Creem un bucle per escriure els resultats amb les dades del vector 'xnums'
    do i = 10000,1000000, 10000
        call montecarlo2(fun2, i, i2, s2)
        write(1,*) i, i2, s2
    end do

    !Apartat 2
    write(1,'(/,A,/)') '#Apartat 2'
    !Creem un bucle per escriure els resultats obtinguts
    do i = 10000, 300000, 10000
        call encerterror(i, a,b,M, fun2, i2, s2)
        write(1,*) i, i2, s2
    end do

end program 

!Subrutina montecarlo
subroutine montecarlo(a, b, fun, n, integral, sigma)
    implicit none
    double precision:: a, b, x, fun, integral, sigma, h, suma, sumaq, mitjq, mitj
    integer:: n, i
    external:: fun
    suma = 0
    sumaq = 0

    !Amb un bucle cridem un nombre aleatori entre 0 i 1
    do i = 1,n
        call random_number(x)
        h = (b - a)*fun((b-a)*x + a) !Fem el canvi de variable per obtenir 'h'
        suma = suma + h !Fem el sumatori
        sumaq = sumaq + h**2    !Fem el sumatori amb h^2 per trobar l'error
    end do

    integral = suma/n   !Calculem la integral amb el valor del sumatori

    !Calculem la mitjana i la mitjana de h^2 per trobar l'error
    mitj = suma/n
    mitjq = sumaq/n

    !Calculem l'error
    sigma = sqrt((mitjq-mitj**2)/n)
    return
end subroutine

!Funció integral 1
double precision function fun1(x)
    implicit none
    double precision:: x

    fun1 = (x**3)*cos(x)**2
    return
end function

!Subrutina acceptrebuig
subroutine acceptrebuig(ndat, xnums, a,b,M, fun)
    implicit none
    double precision:: f, a, b, M, fun, x, x1, x2, p
    double precision:: xnums(ndat), mitjana, suma, var, des
    integer:: i, ndat
    external::fun

    !Càlcul de x aleatòries amb el criteri que ens han donat a la prepràctica
    do i = 1, ndat
        f = 0
        p = 1
        do while (f<p) 
            !Cridem dos nombres aleatoris
            call random_number(x1)
            call random_number(x2)
            !Definim x i p
            x = (b-a)*x1 + a
            p = M*x2
            !Calculem f
            f = fun(x)
        end do 
        !omplim el vector xnums
        xnums(i) = x
    end do

    !Càlcul de la mitjana
    suma = 0
    do i = 1,ndat
        suma = suma + xnums(i)
    end do

    mitjana = suma/ndat

    !Càlcul de la variància
    var = 0
    do i = 1,ndat
        var = var + (xnums(i)-mitjana)**2
    end do
    var = var/(ndat-1)

    !Càlcul de la desviació típica
    des = sqrt(var)

    !Escrivim els resultats al fitxer de dades
    !write(1,*) mitjana, var, des

    return

end subroutine

!Funció de probabilitat per l'apartat b

double precision function p(x)
    implicit none
    double precision:: x, pi, L
    common/dades/pi, L

    p = (1/L)*sin(pi*(x-2*L)/(2*L))**2
    return
end function

!Funció de la integral 2
double precision function fun2(x)
    implicit none
    double precision:: x, pi, L
    common/dades/pi, L

    fun2 = ((1/L)*sin(pi*(x-2*L)/L)**2)
    return
end function

!Subrtutina montecarlo amb els valors del vector xnums
subroutine montecarlo2(fun, n, integral, sigma)
    implicit none
    double precision:: x, fun, integral, sigma, h, suma, sumaq, mitjq, mitj, xnums(1000000)
    integer:: n, i
    external:: fun
    common/vectors/xnums

    suma = 0
    sumaq = 0
    !Creem un bucle per calcular el sumatori
    do i = 1,n
        x = xnums(i) !Definim la x amb els valors del vector 'xnums'
        h = fun(x)  !El valor de 'h' correspon a la imatge de la funció que li hem passat

        !Calculem el sumatori i el sumatori de h^2
        suma = suma + h 
        sumaq = sumaq + h**2
    end do

    !Trobem el valor de la integral a partir del sumatori calculat
    integral = suma/n

    !Calculem la mitjana i la mitjana de h^2
    mitj = suma/n
    mitjq = sumaq/n

    !Trobem l'error
    sigma = sqrt((mitjq-mitj**2)/n)
    return
end subroutine

!Subrutina encerterror
subroutine encerterror(ndat, a,b,M, fun, integral, sigma)
    implicit none
    double precision:: f, a, b, M, fun, x, x1, x2, p
    double precision:: suma, integral, sigma
    integer:: i, ndat
    external::fun

    !Càlcul de x aleatòries amb el criteri que ens han donat a la prepràctica
    do i = 1, ndat
        !Cridem dos nombres aleatoris
        call random_number(x1)
        call random_number(x2)
        !Definim x i p
        x = (b-a)*x1 + a
        p = M*x2
        !Calculem f
        f = fun(x)
        
        !Apliquem el criteri de la teoria
        if (f>p) then
            suma = suma + 1
        end if
    end do

    !Trobem el valor de l'integral i l'error a partir dels punts que es troben a dins la funció
    integral = M*(b-a)*suma/ndat
    sigma = (M*(b-a)/sqrt(ndat*1.d0))*sqrt(suma*1.d0/ndat*(1.d0 - suma*1.d0/ndat))
    return
    
end subroutine