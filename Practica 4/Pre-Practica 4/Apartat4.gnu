set ylabel "error"
set xlabel "h"
set key top left
set xzeroaxis
set yzeroaxis
set logscale y
set format y "10^%L"
set logscale x
set format x "10^%L"
set title "Convergencia entre la massa amb canvi de variables i sense"
set term png
set output "P2-23-24-fig3.png"
plot "P4-23-24-res4.dat" using 1:2 title"Trapezis", "P4-23-24-res4.dat" using 1:3 title"Simpson 3/8","P4-23-24-res4.dat" using 1:4 title"Trapezis canvi de var", "P4-23-24-res4.dat" using 1:5 title"Simpson 3/8 canvi de var"  