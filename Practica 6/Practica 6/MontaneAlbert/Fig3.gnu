set ylabel "I3"
set xlabel "N"

set term png
set output "P6-23-24-fig3.png"

plot  "P6-23-24-res.dat" index 2 using 1:2:3 with errorbars title"I3"