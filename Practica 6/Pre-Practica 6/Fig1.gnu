set ylabel "error"
set xlabel "N"

set term png
set output "P6-23-24-fig1.png"

valor_real1 = pi**3/2
valor_real2 = 6*pi**3 - (905/144)*pi
plot  "P6-23-24-res.dat" index 0 using 1:3 title "Error estimat 1","P6-23-24-res.dat" index 0 using 1:5 title "Error estimat 2","P6-23-24-res.dat" index 0 using 1:($2-valor_real1) title "Error real 1","P6-23-24-res.dat" index 0 using 1:($4-valor_real2) title "Error real 2"