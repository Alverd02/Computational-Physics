set ylabel "x(x_3) [cm]"
set xlabel "x_3 [cm]"
set xzeroaxis
set yzeroaxis
set yrange[0:60]

set term png
set output "P2-23-24-fig2.png"
plot "P2-23-24-res1.dat" using 4:2 title"Pisto 1" with lines, "P2-23-24-res1.dat" using 4:5 title"Pisto 2" with lines