set ylabel "x(t) [cm]"
set xlabel "t [s]"
set xzeroaxis
set yzeroaxis

set term png
set output "P2-23-24-fig1.png"
plot "P2-23-24-res1.dat" using 1:3 title"Pisto 2", "P2-23-24-res1.dat" using 1:4 title"Pisto 3"