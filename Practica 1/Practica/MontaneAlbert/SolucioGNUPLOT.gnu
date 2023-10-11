#------------------------------------------------------------------------------------------
#--    PLANTILLA d'script GNUPLOT per fer grafiques més traballades
#--    Les opcions no activades estan precedides d'un # a fi i efecte de no oblidar-les.
#--
#--    R. Mayol (9/2018)
#------------------------------------------------------------------------------------------

set zero 1.e-20               # Agafa com zero quantitats inferiors
#set encoding  iso_8859_15     # Permet l'ús d'accents als  MAS 
#set encoding  iso_8859_1      # Permet l'ús d'accents als PC's    
#set encoding  CP1250          # Permet l'ús d'accents als MS-Windows
#set locale "es_ES.utf8"

#-------------------------------------------------------------------
#------ A vegades la grafica no acaba de encabir-se a la finestra
#------ amb això es controla una mica més el tamany 
#-------------------------------------------------------------------

set lmargin at screen 0.15 
set rmargin at screen 0.9
set bmargin at screen 0.15
set tmargin at screen 0.8 
 
#-------------------------------------
#------ Legenda 
#-------------------------------------

#set key left                    # Llegenda al costat esquerre 
#set key font "Arial, 20"        # Font de la legenda

#set key spacing 1.5 R font  "Arial, 18" # Espaiat de línies i tamany font a la legenda

set key at 200, 1.e4             # Llegenda en aquestes coordenades
#set key box lt -1 lw 1           # Recepte per posar una caixa a la legenda
#set key spacing 1.5 R            # Espaiat de línies

#--------------------------------------
#-- X-Axis
#--------------------------------------

set xlabel "N (Línies) "      # Etiqueta de l'eix X
set xrange [0.:220.]          # Limits de l'eix X

#set xtics 0,4,40              
#set mxtics 5                 # Numero de subtics per cada xtic
#set xtics border mirror 10

set grid xtics                 # Linies de grid X fix
#set xtics  font "Arial, 20"    # Font per les etiquetes situades als xtics
#set xlabel font "Arial, 20"     
#set logscale x                 # Activa l'escala logaritmica

#--------------------------------------
#-- Y-Axis
#--------------------------------------


#set ytics  font "Arial, 20"
#set ylabel font "Arial, 20"
#set yrange [5.e5:1.e6]
#set ytics 5.e5,1.e5,1.e6
#set mytics 5
#set ytics border mirror 1

set ylabel "S^{3}_{N} (porcions) "
set grid ytics
set format y "10^%L"        # Les etiquetes logaritmiques són millors.
set logscale y

#--------------------------------------
#-- Labels and Title
#--------------------------------------

#set label 1 "bla bla " at screen 0.5,0.69 center font "Helvetica, 24"

set title 'P1-1819-res1'
set title font "Arial, 16"

#set label "S^{3}_{N} = Número de porcions amb N-2 pizzes" at screen 0.35, 0.83 font "Arial,18"

#--------------------------------
#--- PLOT
#--------------------------------

plot "P1-1819-res1.dat" using 1:2            with lines linewidth 2 linecolor 1 t'S^{3}_{N}' , \
     "P1-1819-res1.dat" using 1:(($1**3)*5./9.)  with lines linewidth 2 linecolor 2 t'S^{asim}_{N}'




pause -1    # Espera fins que es toqui una tecla


#------ Sortida a arxiu png 
set terminal png enhanced

set output "./P4-1819-fig1.png"
replot

