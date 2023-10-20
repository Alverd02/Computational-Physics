set ylabel "x(t) [cm]"
set xlabel "t [s]"
set xzeroaxis
set yzeroaxis
set yrange[0:45]

set term png
set output "P2-23-24-fig1-c.png"
plot "P2-23-24-res1-c.dat" using 1:2 title"Pisto 1", "P2-23-24-res1-c.dat" using 1:4 title"Pisto 3","P2-23-24-res1-c.dat" using 1:6 title"Pisto 5"