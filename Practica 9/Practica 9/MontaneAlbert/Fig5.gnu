set term png
set output "P9-23-24-fig5.png"
set view map
set size ratio -1
set tics out
set tics nomirror
set yrange [0:45.5]
set xrange [0:33.5]
set xlabel "y (cm)"
set ylabel "x (cm)"
set pm3d map

splot "res.dat" index 9 using 2:1:3 with pm3d notitle