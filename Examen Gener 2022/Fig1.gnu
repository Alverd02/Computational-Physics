set ylabel "p(x)"
set xlabel "x "


set term png
set output "figE1.png"

p(x) = exp(-x**2/(2*(sqrt(2))**2))/(sqrt(2)*sqrt(2*pi))
plot  p(x),"apartat1a.dat" index 0 using 1:2:3 with errorbars title"histograma"