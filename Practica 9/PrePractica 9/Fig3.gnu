set term png
set output "figura2.png"
set view map
set size ratio -1
set tics out
set tics nomirror
set yrange [0:18.5]
set xrange [0:31]
set xlabel "y (cm)"
set ylabel "x (cm)"
set pm3d map

splot "mapa.dat" using 2:1:3 with pm3d notitle










