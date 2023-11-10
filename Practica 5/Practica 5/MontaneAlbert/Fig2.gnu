set ylabel "g(x)  [micrometres^-1]"
set xlabel "x [micrometres]"


set term png
set output "P2-23-24-fig2.png"

f(x) = (exp(-(x**2)/(2*3**2))/sqrt(2*pi*3**2))
plot  f(x),"P5-23-24-res.dat" index 2 using 1:2:3 with errorbars title"histograma"