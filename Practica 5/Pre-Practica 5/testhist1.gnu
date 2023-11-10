set term png
set output "P2-23-24-fig1.png"

f(x) = (125*exp(pi)*x**2*(sin(x))**2*exp(-abs(x)))/(4*(68*exp(pi)-70*pi-25*pi**2-68))
plot  f(x),"P5-23-24-res.dat" index 1 using 1:2:3 with errorbars title"histograma"