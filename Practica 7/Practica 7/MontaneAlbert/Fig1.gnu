set ylabel "theta [rad]"
set xlabel "t [s]"
set key bmargin center horizontal box

set term png
set output "P6-23-24-fig1.png"

f(x) = 0.025*cos(x*sqrt(10.44/1.07))
plot  "P7-23-24-res.dat" index 0 using 1:2 title "Euler" with lines, f(x) title "f(x) aprox","P7-23-24-res.dat" index 0 using 1:3 with lines title "Euler segon ordre" 