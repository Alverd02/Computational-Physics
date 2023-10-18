set ylabel "x(t) [cm]"
set xlabel "t [s]"
set xzeroaxis
set yzeroaxis
set yrange[0:60]


set term png
set output "P2-23-24-fig3.png"
plot "P2-23-24-res2.dat" title"Pisto 2 interpolat" with lines, "P2-23-24-res1.dat" using 1:3 title"Pisto 2 sense interpolació"