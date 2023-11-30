program p7
    implicit none
    double precision:: l, g, m, y0, p0, t0, tf, tn, pi, osci, h, t, ekin, epot, etot
    double precision,dimension(:), allocatable :: y, p
    integer:: k, i
    external:: osci
    common/dades/l,g,m, pi

    !Definim les dades que utilitzarem al llarg de la pràctica
    m = 0.98
    l = 1.07
    g = 10.44
    pi = 4*atan(1.d0)
    tn = 2*pi*sqrt(l/g)

    !Apartat a - petites oscilacions

    !En aquest apartat i els dos següents farem 1500 passos
    k = 1500

    !Definim la dimensió dels vectors y i p, que són la posició i la velocitat respectivament
    allocate(y(k), p(k))
    
    !Definim les condicions inicals, l'interval de temps i el salt
    y0 = 0.025
    p0 = 0
    t0 = 0
    tf =7*tn
    h = (tf - t0)/(k*1.d0)

    open(1, file= 'P7-22-23-res.dat')

    write(1,'(/,A,/)') '#Apartat a'
    write(1,'(/,A,/)') '#temps, phi i dphi obtinguts amb Euler'

    !Cridem la subrutina 'euler' per trobar els valors de y i de p
    call euler(y0, p0, t0, tf, k, osci, y, p)

    !Escrivim els temps, les posicions i les velocitats obtingudes al fitxer de dades
    do i = 1, k
        t = t0 + h*i
        write(1,'(3F25.14)') t, y(i), p(i)
    end do

    write(1,'(/,A,/)') '#temps, phi i dphi obtinguts amb RK2'

    !Cridem la subrutina rk2
    call rk2(y0,p0,t0,tf,k,osci,y,p)

    !Escrivim els resultats al fitxer de dades
    do i = 1, k
        t = t0 + h*i
        write(1,'(3F25.14)') t, y(i), p(i)
    end do

    !Apartat b - oscil·lacions grans
    y0 = pi - 0.15
    p0 = 0

    write(1,'(/,A,/)') '#Apartat b'
    write(1,'(/,A,/)') '#temps, phi i dphi obtinguts amb Euler'

    !Cridem la subrutina 'euler' per trobar els valors de y i de p
    call euler(y0, p0, t0, tf, k, osci, y, p)

    !Escrivim els temps, les posicions i les velocitats obtingudes al fitxer de dades
    do i = 1, k
        t = t0 + h*i
        write(1,'(3F25.14)') t, y(i), p(i)
    end do

    write(1,'(/,A,/)') '#temps, phi i dphi obtinguts amb RK2'

    !Cridem la subrutina
    call rk2(y0,p0,t0,tf,k,osci,y,p)

    !Escrivim els temps, les posicions i les velocitats obtingudes al fitxer de dades
    do i = 1, k
        t = t0 + h*i
        write(1,'(3F25.14)') t, y(i), p(i)
    end do

    !Apartat c
    write(1,'(/,A,/)') '#Apartat c'
    
    !Definim les condicions inicials
    y0 = pi - 0.025
    p0 = 0.12

    !Cridem la subrutina 'euler' per trobar els valors de y i de p
    call euler(y0, p0, t0, tf, k, osci, y, p)

    write(1,'(/,A,/)') '#temps, energia cinètica i energia potencial calculades amb Euler'

    !amb l'ajuda de les funcions 'ekin' i 'epot' calculem els valors de l'energia mecància
    !Escrivim els resultats al fitxer
    do i = 1, k
        t = t0 + h*i
        etot= ekin(p(i)) + epot(y(i))
        write(1,'(4F25.14)') t, ekin(p(i)), epot(y(i)), etot
    end do

    !Cridem la subrutina 'rk2' per trobar els valors de y i de p  
    call rk2(y0,p0,t0,tf,k,osci,y,p)

    !Calculem els valors de l'energia mecànica i escrivim els resultats al fitxer de dades
    write(1,'(/,A,/)') '#temps, energia cinètica i energia potencial calculades amb Adams-Bashforth'
    do i = 1, k
        t = t0 + h*i
        etot= ekin(p(i)) + epot(y(i))
        write(1,'(4F25.14)') t, ekin(p(i)), epot(y(i)), etot
    end do

    !Canviem la dimensió dels vectors
    deallocate(y,p)

    !Apartat d
    !En aquest apartat fem servir 6000 passos
    k = 6000

    !Redefinim la dimensió dels vectors
    allocate(y(k), p(k))

    !Definim les condicions inicials, l'interval de temps i el salt
    y0 = 0
    p0 = 2*sqrt(g/l) + 0.05
    tf = 15*tn
    h = (tf-t0)/(k*1.d0)

    write(1,'(/,A,/)') '#Apartat d - primer cas'
    write(1,'(/,A,/)') '#temps, phi i dphi obtinguts amb RK2'

    !Cridem la subrutina 'rk2' per calcular les posicions i les velocitats
    call rk2(y0, p0, t0, tf, k, osci, y, p)

    !Escrivim els resultats al fitxer
    do i = 1, k
        t = t0 + h*i
        write(1,'(3F25.14)') t, y(i), p(i)
    end do

    write(1,'(/,A,/)') '#Apartat d - segon cas'
    write(1,'(/,A,/)') '#temps, phi i dphi obtinguts amb RK2'

    !Definim la velocitat inicial (segona opció)
    p0 = 2*sqrt(g/l) - 0.05

    !Cridem la subrutina 'rk2' per calcular les posicions i les velocitats
    call rk2(y0, p0, t0, tf, k, osci, y, p)

    !Escrivim els resultats al fitxer
    do i = 1, k
        t = t0 + h*i
        write(1,'(3F25.14)') t, y(i), p(i)
    end do

    deallocate(y,p)

    !apartat e
    write(1,'(/,A,/)') '#Apartat 3'

    !Definim les condicions inicials i l'interval de temps
    y0=2.87
    p0 = 0
    t0 = 0
    tf = 11*tn

    !300 passos
    k=300

    !Amb els passos calculem el salt i redefinim les dimensions dels vectors
    h = (tf-t0)/(k*1.d0)
    allocate(y(k), p(k))

    !Cridem la subrutina 'rk2'
    call rk2(y0,p0,t0,tf,k,osci,y,p)

    write(1,'(/,A,/)') '#temps, energia cinètica i energia potencial calculades amb 300 passos'

    !Calculem l'energia i escrivim els resultats al fitxer
    do i = 1, k
        t = t0 + h*i
        etot= ekin(p(i)) + epot(y(i))
        write(1,'(2F25.14)') t, etot
    end do

    deallocate(y,p)

    !Repetim el procés per 550, 1000 i 20000 passos
    !550 passos
    k=550
    h = (tf-t0)/(k*1.d0)
    allocate(y(k), p(k))

    call rk2(y0,p0,t0,tf,k,osci,y,p)

    write(1,'(/,A,/)') '#temps, energia cinètica i energia potencial calculades amb 650 passos'
    do i = 1, k
        t = t0 + h*i
        etot= ekin(p(i)) + epot(y(i))
        write(1,'(2F25.14)') t, etot
    end do

    deallocate(y,p)

    !1000 passos
    k=1000
    h = (tf-t0)/(k*1.d0)
    allocate(y(k), p(k))

    call rk2(y0,p0,t0,tf,k,osci,y,p)

    write(1,'(/,A,/)') '#temps, energia cinètica i energia potencial calculades amb 1250 passos'
    do i = 1, k
        t = t0 + h*i
        etot= ekin(p(i)) + epot(y(i))
        write(1,'(2F25.14)') t, etot
    end do

    deallocate(y,p)

    !20000 passos
    k=20000
    h = (tf-t0)/(k*1.d0)
    allocate(y(k), p(k))

    call rk2(y0,p0,t0,tf,k,osci,y,p)

    write(1,'(/,A,/)') '#temps, energia cinètica i energia potencial calculades amb 25000 passos'
    do i = 1, k
        t = t0 + h*i
        etot= ekin(p(i)) + epot(y(i))
        write(1,'(2F25.14)') t, etot
    end do

    deallocate(y,p)

    close(1)
    
end program 

!Subrutina euler
subroutine euler(y0, p0, t0, tf, k, fun, y, p)
    implicit none
    double precision:: y0, p0, t0, tf, fun, y(k), p(k), h
    integer:: k, i
    external:: fun

    !Calculem el salt i definim el primer valor dels vectors posició i velocitat, que corresponen als de les condicions inicials
    h = (tf-t0)/(k*1.d0)
    y(1) = y0
    p(1) = p0

    !Calculem els valors de les següents posicions i velocitats
    do i = 2, k
        y(i) = y(i-1) + h*p(i-1)
        p(i) = p(i-1) + h*fun(y(i-1))
    end do
    return

end subroutine

subroutine rkv(x0, dim, yyin, h, yyout, tipus)
    implicit none
    integer:: tipus, dim
    double precision:: x0, h
    double precision:: yyin(dim), yyout(dim), k1(dim), k2(dim), k3(dim), k4(dim)

    !Per trobar el valor de k1, k2, k3 i k4 cridem la subrutina que ens retorna 
    !la imatge de la funció que ens proporciona l'enunciat. Aquesta funció és un 
    !vector de dimensió 'dim'

    call fun(dim, x0, yyin, k1)
    call fun(dim, x0 + h/2.d0,  yyin + (h/2.d0)*k1, k2)
    call fun(dim, x0 + h/2.d0,  yyin + (h/2.d0)*k2, k3)
    call fun(dim, x0 + h, yyin + h*k3, k4)

    !Segons el tipus de mètode que vulguem (RK2, RK3 o RK4), trobarem els valors
    !Del següent punt a integrar
    
    if (tipus.eq.2) then
        yyout = yyin + h*k2
        return

    else if (tipus.eq.3) then
        yyout = yyin + (h/6.d0)*(k1 + 4*k2 + k3)
        return

    elseif (tipus.eq.4) then
        yyout = yyin + (h/6.d0)*(k1 + 2*k2 + 2*k3 + k4)
        return
    endif

end subroutine

!Funció 'osci', equació del moviment del pèndol simple
double precision function osci(x)
    implicit none
    double precision:: x, l, g, m, pi
    common/dades/l, g, m, pi

    osci = (-g/l)*dsin(x)

    return
end function

!Funció 'ekin', equació que ens dona l'energia cinètica
double precision function ekin(x)
    implicit none
    double precision:: x, l, g, m, pi
    common/dades/l, g, m, pi

    ekin = (1.d0/2)*m*(x**2)*(l**2)

    return
end function

!Funció 'epot', equació que ens dona l'energia potencial
double precision function epot(x)
    implicit none
    double precision:: x,l, g, m, pi
    common/dades/l, g, m, pi

    epot = -m*g*l*dcos(x)

    return
end function