set term png
set output "P3-23-24-fig3.png"

set multiplot layout 1,2

set ylabel "f(x)'"
set xlabel "x"
set xzeroaxis
set yzeroaxis
set xrange[0:2*pi]
set title "25 punts"


plot "P3-23-24-res3-n25.dat" using 1:3 title"Derivada eproximada","P3-23-24-res3-n25.dat" using 1:4 title"Derivada exacte"

set ylabel "f(x)'"
set xlabel "x"
set xzeroaxis
set yzeroaxis
set xrange[0:2*pi]
set title "230 punts"


plot "P3-23-24-res3-n230.dat" using 1:3 title"Derivada eproximada","P3-23-24-res3-n230.dat" using 1:4 title"Derivada exacte"

unset multiplot