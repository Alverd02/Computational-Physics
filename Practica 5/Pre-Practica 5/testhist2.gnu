set term png
set output "P2-23-24-fig2.png"

f(x) = (6/7)*exp(-(6/7)*x)
plot  f(x),"P5-23-24-res.dat" index 3 using 1:2:3 with errorbars title"histograma"