set ylabel "Número de quarks"
set xlabel "N"

set term png
set output "P6-23-24-fig1.png"

plot  "P6-23-24-res.dat" index 0 using 1:2:3 with errorbars title"n_u","P6-23-24-res.dat" index 0 using 1:4:5 with errorbars title"n_d"