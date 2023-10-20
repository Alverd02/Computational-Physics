set ylabel "x(x_1) [cm]"
set xlabel "log(x_1)"
set xzeroaxis
set yzeroaxis
set yrange[10:30]
set logscale x
set format x "10^%L"
 
set term png
set output "P2-23-24-fig2-c.png"
plot "P2-23-24-res1-c.dat" using 2:4 title"Pisto 3" with lines, "P2-23-24-res1-c.dat" using 2:6 title"Pisto 5" with lines