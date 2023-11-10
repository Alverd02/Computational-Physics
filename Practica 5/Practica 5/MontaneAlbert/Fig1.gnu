set ylabel "p(x) [nm^-1]"
set xlabel "x [nm]"
set yrange [0:0.1]

set term png
set output "P2-23-24-fig1.png"

f(x) = sin(x/4)**2/(4*pi)
plot  f(x),"P5-23-24-res.dat" index 0 using 1:2:3 with errorbars title"histograma"