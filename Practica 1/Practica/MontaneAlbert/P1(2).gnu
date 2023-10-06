set xlabel "N"
set ylabel "S^8_N/S^(asim)_N"
set xzeroaxis
set yzeroaxis

set term png
set output "P1-23-24-fig2.png"

plot "P1-23-24-res2.dat"