set ylabel "I2"
set xlabel "N"

set term png
set output "P6-23-24-fig2.png"

plot  "P6-23-24-res.dat" index 1 using 1:2:3 with errorbars title"I2"