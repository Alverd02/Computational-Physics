set ylabel "x(t) [cm]"
set xlabel "t [s]"
set xzeroaxis
set yzeroaxis
set yrange[10:30]


set term png
set output "P2-23-24-fig3-c.png"
plot "P2-23-24-res2-c.dat" using 1:3 title"Pisto 4 interpolat lineal" with lines, "P2-23-24-res1-c.dat" using 1:5 title"Pisto 4","P2-23-24-res2-c.dat" using 1:2 title"Pisto 4 interpolat ordre 0" 