set ylabel "f(x)'"
set xlabel "x"
set xzeroaxis
set yzeroaxis
set xrange[0:2*pi]



set term png
set output "P3-23-24-fig3.png"

plot "P2-23-24-res3-n25.dat" using 1:3 title "Pisto 1"