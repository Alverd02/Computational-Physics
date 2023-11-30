       PROGRAM PREPRACTICA7
       IMPLICIT NONE
       EXTERNAL PHIPP,PHIPPAPROX

       INTEGER I,N,J,PAS(4)
       DOUBLE PRECISION M,G,L,T,TN,WN,PI,PHI,DPHI,PHI0,DPHI0,H,
     1 PHIM1,DPHIM1,PHIA,DPHIA,PHI0A,DPHI0A,PHIM1A,DPHIM1A,
     1 EPOT,ECIN
       COMMON/DADES/M,L,G

       M=0.51D0
       L=45.D-2
       PI=4*DATAN(1.D0)
       G=8.87D0
       WN=DSQRT(G/L)
       TN=2*PI/WN

C a) Oscil·lacions grans:
       OPEN(14,FILE="P7_1920_res.dat")
C mètode d'Euler
       WRITE(14,*) "t,    Phi,   dPhi/dt    (mètode d'Euler)"
       N=1800
       H=5*TN/N
       PHI0=PI-0.02D0
       DPHI0=0.D0
       WRITE(14,15) 0.D0,PHI0,DPHI0
15     FORMAT(3(E14.6,2X))
       DO I=1,N,1
        T=H*I
C calcula cada valor nou agafant l'anterior com a inicial
        CALL EULER(PHI,DPHI,PHI0,DPHI0,PHIPP,H)
        WRITE(14,15) T,PHI,DPHI
        PHI0=PHI
        DPHI0=DPHI
       ENDDO


C mètode d'Euler millorat
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "t,    Phi,   dPhi/dt    (Euler millorat)"
C condicions inicials (en calen 4 en aquest cas): 
       PHIM1=PI-0.02D0
       DPHIM1=0.D0
       WRITE(14,15) 0.D0,PHIM1,DPHIM1
C el segon valor es calcula amb el mètode d'Euler
       CALL EULER(PHI0,DPHI0,PHIM1,DPHIM1,PHIPP,H)
       WRITE(14,15) H,PHI0,DPHI0
       DO I=2,N,1
        T=H*I
C s'agafa com a inicial el valor de dos iteracions enrere
        CALL EULERM(PHI,DPHI,PHIM1,DPHIM1,PHIPP,H)
        WRITE(14,15) T,PHI,DPHI
        PHIM1=PHI0
        DPHIM1=DPHI0
        PHI0=PHI
        DPHI0=DPHI
       ENDDO


C b) Oscil·lacions petites:

C mètode d'Euler

       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "t, Phi, Phi', Phi aprox, Phi' aprox (Euler)"
       N=1300
       H=5*TN/N
       PHI0=0.04D0
       DPHI0=0.D0
       PHI0A=0.04D0
       DPHI0A=0.D0
       WRITE(14,16) 0.D0,PHI0,DPHI0,PHI0A,DPHI0A
16     FORMAT(5(E11.4,2X))
       DO I=1,N,1
        T=H*I
C calcula cada valor nou agafant l'anterior com a inicial
        CALL EULER(PHIA,DPHIA,PHI0A,DPHI0A,PHIPPAPROX,H)
        CALL EULER(PHI,DPHI,PHI0,DPHI0,PHIPP,H)
        WRITE(14,16) T,PHI,DPHI,PHIA,DPHIA
        PHI0=PHI
        DPHI0=DPHI
        PHI0A=PHIA
        DPHI0A=DPHIA
       ENDDO


C mètode d'Euler millorat

       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "t, Phi, Phi', Phi aprox, Phi' aprox (millorat)"
C condicions inicials (en calen 4 en aquest cas): 
       PHIM1=0.04D0
       PHIM1A=0.04D0
       DPHIM1=0.D0
       DPHIM1A=0.D0
       WRITE(14,16) 0.D0,PHIM1,DPHIM1,PHIM1A,DPHIM1A
C el segon valor es calcula amb el mètode d'Euler
       CALL EULER(PHI0,DPHI0,PHIM1,DPHIM1,PHIPP,H)
       CALL EULER(PHI0A,DPHI0A,PHIM1A,DPHIM1A,PHIPPAPROX,H)
       WRITE(14,16) H,PHI0,DPHI0,PHI0A,DPHI0A
       DO I=2,N,1
        T=H*I
C s'agafa com a inicial el valor de dos iteracions enrere
        CALL EULERM(PHI,DPHI,PHIM1,DPHIM1,PHIPP,H)
        CALL EULERM(PHIA,DPHIA,PHIM1A,DPHIM1A,PHIPPAPROX,H)
        WRITE(14,16) T,PHI,DPHI,PHIA,DPHIA
        PHIM1=PHI0
        DPHIM1=DPHI0
        PHIM1A=PHI0A
        DPHIM1A=DPHI0A
        PHI0=PHI
        DPHI0=DPHI
        PHI0A=PHIA
        DPHI0A=DPHIA
       ENDDO


C c) Energia------------------------------------------------------


C mètode d'Euler
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "t,    K(t),   V(t),   E(t)    (Phi0=1, Euler)"
       N=2500
       H=5*TN/N
       PHI0=1.D0
       DPHI0=0.D0
       WRITE(14,17) 0.D0,ECIN(DPHI0),EPOT(PHI0),ECIN(DPHI0)+EPOT(PHI0)
17     FORMAT(4(E14.6,2X))
       DO I=1,N,1
        T=H*I
C calcula cada valor nou agafant l'anterior com a inicial
        CALL EULER(PHI,DPHI,PHI0,DPHI0,PHIPP,H)
        WRITE(14,17) T,ECIN(DPHI),EPOT(PHI),ECIN(DPHI)+EPOT(PHI)
        PHI0=PHI
        DPHI0=DPHI
       ENDDO


C mètode d'Euler millorat
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "t,    K(t),   V(t),    E(t)    (Phi0=1, millorat)"
C condicions inicials (en calen 4 en aquest cas): 
       PHIM1=1.D0
       DPHIM1=0.D0
       WRITE(14,17) 0.D0,ECIN(PHIM1),EPOT(DPHIM1),
     1 ECIN(PHIM1)+EPOT(DPHIM1)
C el segon valor es calcula amb el mètode d'Euler
       CALL EULER(PHI0,DPHI0,PHIM1,DPHIM1,PHIPP,H)
       WRITE(14,17) H,ECIN(PHI0),EPOT(DPHI0),
     1 ECIN(PHI0)+EPOT(DPHI0)
       DO I=2,N,1
        T=H*I
C s'agafa com a inicial el valor de dos iteracions enrere
        CALL EULERM(PHI,DPHI,PHIM1,DPHIM1,PHIPP,H)
        WRITE(14,17) T,ECIN(DPHI),EPOT(PHI),ECIN(DPHI)+EPOT(PHI)
        PHIM1=PHI0
        DPHIM1=DPHI0
        PHI0=PHI
        DPHI0=DPHI
       ENDDO

C canviem les condicions inicials i repetim el mateix procedimant

C mètode d'Euler
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "t,    K(t),   V(t),   E(t)    (Phi0=PI-0,02, Euler)"
       N=2500
       H=5*TN/N
       PHI0=PI-0.02D0
       DPHI0=0.D0
       WRITE(14,17) 0.D0,ECIN(DPHI0),EPOT(PHI0),ECIN(DPHI0)+EPOT(PHI0)
       DO I=1,N,1
        T=H*I
C calcula cada valor nou agafant l'anterior com a inicial
        CALL EULER(PHI,DPHI,PHI0,DPHI0,PHIPP,H)
        WRITE(14,17) T,ECIN(DPHI),EPOT(PHI),ECIN(DPHI)+EPOT(PHI)
        PHI0=PHI
        DPHI0=DPHI
       ENDDO


C mètode d'Euler millorat
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "t,   K(t),  V(t),  E(t) (Phi0=PI-0,02, millorat)"
C condicions inicials (en calen 4 en aquest cas): 
       PHIM1=PI-0.02D0
       DPHIM1=0.D0
       WRITE(14,17) 0.D0,ECIN(PHIM1),EPOT(DPHIM1),
     1 ECIN(PHIM1)+EPOT(DPHIM1)
C el segon valor es calcula amb el mètode d'Euler
       CALL EULER(PHI0,DPHI0,PHIM1,DPHIM1,PHIPP,H)
       WRITE(14,17) H,ECIN(PHI0),EPOT(DPHI0),
     1 ECIN(PHI0)+EPOT(DPHI0)
       DO I=2,N,1
        T=H*I
C s'agafa com a inicial el valor de dos iteracions enrere
        CALL EULERM(PHI,DPHI,PHIM1,DPHIM1,PHIPP,H)
        WRITE(14,17) T,ECIN(DPHI),EPOT(PHI),ECIN(DPHI)+EPOT(PHI)
        PHIM1=PHI0
        DPHIM1=DPHI0
        PHI0=PHI
        DPHI0=DPHI
       ENDDO


C d) Transició

C faig servir les mateixes variables que en l'aproximació per petites
C oscil·lacions, on ara PHI és per la condició inicial amb el + i 
C PHIA amb el -.
       N=2100
       H=11*TN/N
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "t, Phi+, Phi'+, Phi-, Phi'-"
C condicions inicials: 
       PHIM1=0.D0
       PHIM1A=0.D0
       DPHIM1=2*DSQRT(G/L)+0.03D0
       DPHIM1A=2*DSQRT(G/L)-0.03D0
       WRITE(14,16) 0.D0,PHIM1,DPHIM1,PHIM1A,DPHIM1A
C el segon valor es calcula amb el mètode d'Euler
       CALL EULER(PHI0,DPHI0,PHIM1,DPHIM1,PHIPP,H)
       CALL EULER(PHI0A,DPHI0A,PHIM1A,DPHIM1A,PHIPP,H)
       WRITE(14,16) H,PHI0,DPHI0,PHI0A,DPHI0A
       DO I=2,N,1
        T=H*I
C s'agafa com a inicial el valor de dos iteracions enrere
        CALL EULERM(PHI,DPHI,PHIM1,DPHIM1,PHIPP,H)
        CALL EULERM(PHIA,DPHIA,PHIM1A,DPHIM1A,PHIPP,H)
        WRITE(14,16) T,PHI,DPHI,PHIA,DPHIA
        PHIM1=PHI0
        DPHIM1=DPHI0
        PHIM1A=PHI0A
        DPHIM1A=DPHI0A
        PHI0=PHI
        DPHI0=DPHI
        PHI0A=PHIA
        DPHI0A=DPHIA
       ENDDO


C d) Convergència del mètode
       PAS=[600,1300,2600,15000]
       WRITE(14,*)
       WRITE(14,*)
       WRITE(14,*) "t, Phi+, Phi'+, Phi-, Phi'-"
       DO J=1,4,1     
        N=PAS(J)
        H=12*TN/N
C condicions inicials (en calen 4 en aquest cas): 
        PHIM1=2.8D0
        DPHIM1=0.D0
        WRITE(14,18) 0.D0,ECIN(PHIM1)+EPOT(DPHIM1)
C el segon valor es calcula amb el mètode d'Euler
       CALL EULER(PHI0,DPHI0,PHIM1,DPHIM1,PHIPP,H)
       WRITE(14,18) H,ECIN(PHI0)+EPOT(DPHI0)
       DO I=2,N,1
        T=H*I
C s'agafa com a inicial el valor de dos iteracions enrere
        CALL EULERM(PHI,DPHI,PHIM1,DPHIM1,PHIPP,H)
        WRITE(14,18) T,ECIN(DPHI)+EPOT(PHI)
        PHIM1=PHI0
        DPHIM1=DPHI0
        PHI0=PHI
        DPHI0=DPHI
       ENDDO
       WRITE(14,*)
       WRITE(14,*)
       ENDDO
18     FORMAT(2(E14.6,2X))







       CLOSE(14)
       CALL SYSTEM("gnuplot fig.gnu")

       END PROGRAM



C DEFINICIÓ DE LES SUBRUTINES----------------------------
C Mètode d'Euler
       SUBROUTINE EULER(Y,DY,Y0,DY0,FUN,H)
C mètode d'Euler per calcular el valors de Y(X) i DY(X)=Y'(X) amb un pas de 
C H partint de la segona derivada FUN=Y''(X) i les condicions inicials Y0,DY0.
       IMPLICIT NONE
       DOUBLE PRECISION Y,Y0,H,DY0,DY,DDY0
       CALL FUN(Y0,DDY0)
       DY=DY0+H*DDY0
       Y=Y0+H*DY0
       RETURN
       END

C mètode d'Euler millorat
       SUBROUTINE EULERM(Y,DY,YM1,DYM1,FUN,H)
C mètode d'Euler millorat per calcular el valors de Y(X) i DY(X)=Y'(X) amb un pas de 
C H partint de la segona derivada FUN=Y''(X) i les condicions inicials YM1,DYM1, que en 
C aquest cas corresponen al terme n-1 per calcular el terme n+1.
       IMPLICIT NONE
       DOUBLE PRECISION Y,YM1,H,DYM1,DY,DDYM1
       CALL FUN(YM1,DDYM1)
       DY=DYM1+2*H*DDYM1
       Y=YM1+2*H*DYM1
       RETURN
       END


C Equació diferencial del moviment (derivada segona)
       SUBROUTINE PHIPP(PHI,DDPHI)
C retorna la segona derivada de PHI(t) al punt PHI1
       IMPLICIT NONE
       DOUBLE PRECISION M,G,L,PHI,DDPHI
       COMMON/DADES/M,L,G
       DDPHI=-G*DSIN(PHI)/L
       RETURN
       END

C Equació aproximada per oscil·lacions petites
       SUBROUTINE PHIPPAPROX(PHI,DDPHIAP)
C retorna la segona derivada de PHI(t) al punt PHI1
       IMPLICIT NONE
       DOUBLE PRECISION M,G,L,PHI,DDPHIAP
       COMMON/DADES/M,L,G
       DDPHIAP=-G*(PHI)/L
       RETURN
       END

C DEFINICIÓ DE LES FUNCIONS
C Energia cinètica
       DOUBLE PRECISION FUNCTION ECIN(DPHI)
       IMPLICIT NONE
       DOUBLE PRECISION M,L,G,DPHI
       COMMON/DADES/M,L,G
       ECIN=0.5D0*M*DPHI**2*L**2
       RETURN
       END

C Energia potencial
       DOUBLE PRECISION FUNCTION EPOT(PHI)
       IMPLICIT NONE
       DOUBLE PRECISION M,L,G,PHI
       COMMON/DADES/M,L,G
       EPOT=-M*G*L*DCOS(PHI)
       RETURN
       END
