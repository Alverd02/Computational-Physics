set term png
set output "P1-23-24-fig2.png"
set xrange [0:20]
plot 8-6/x+6/x**2
plot "P1-23-24(2).dat" 