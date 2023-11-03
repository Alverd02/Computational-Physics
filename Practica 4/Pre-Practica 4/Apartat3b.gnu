set ylabel "error"
set xlabel "h"
set xzeroaxis
set key top left
set yzeroaxis
set logscale y
set format y "10^%L"
set logscale x
set format x "10^%L"
set title "Convergencia apartat b"
set term png
set output "P2-23-24-fig2.png"
plot "P4-23-24-res3.dat" using 1:2 title"Trapezis teoric" with lines, "P4-23-24-res3.dat" using 1:3 title"Simpson teoric 3/8" with lines,"P4-23-24-res3.dat" using 1:4 title"Trapezis", "P4-23-24-res3.dat" using 1:5 title"Simpson 3/8"  