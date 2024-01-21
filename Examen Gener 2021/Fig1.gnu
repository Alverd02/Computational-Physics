set ylabel "T [s]"
set xlabel "Theta0 [rad] "

set key bmargin center horizontal box

set term png
set output "figE1.png"

plot  "apartat1a.dat" index 0 using 1:2 title"Metode de Trapezis" pointtype 5,"apartat1a.dat" index 0 using 1:3:4 with errorbars title"Metode de MonteCarloCru"