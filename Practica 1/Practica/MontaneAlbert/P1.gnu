set xlabel "N"
set ylabel "S^8_N"
set xzeroaxis
set yzeroaxis

set term png
set output "P1-23-24-fig1.png"
set multiplot layout 1,2

plot "P1-23-24-res1.dat", x**3/5 
plot "P1-23-24-res1.dat" using 1:log(2), log(x**3/5) 
unset multiplot