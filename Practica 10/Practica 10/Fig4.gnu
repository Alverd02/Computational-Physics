set term gif size 1000,600 animate  delay 10 loop 0 optimize
set output "P10-1819-fig4.gif"

# leeremos el numero de bloques de manera automatica
datafile ="apartat2c"
stats datafile
numblo=STATS_blocks

# pone el titulo
set title "T(x,t)"

# fijamos los rangos del canvas, para que no se reajusten
set xrange[0:1.5]
set yrange[0:280]

#fijamos los títulos
set xlabel "x [m]"
set ylabel "T(x,t) [ºC]"

#bucle

do for[i=0:numblo-2:1]{
#escribe la etiqueta del tiempo
set label 2 sprintf('Time: %9.3f (sec)',i*0.002) at 10,160 right front font 'Verdana,12'

#dibuja
plot datafile index i u 2:3 w l lw 4 t"T(x,t)","" index 0 u 2:3 t"T_0(x)" w l lw 3

#borra la etiqueta
unset label 2
}
