set ylabel "error"
set xlabel "N"

set term png
set output "P6-23-24-fig1.png"

plot  "P6-23-24-res.dat" index 0 using 1:3 title "Error estimat 1","P6-23-24-res.dat" index 0 using 1:6 title "Error estimat 2","P6-23-24-res.dat" index 0 using 1:4 title "Error real 1","P6-23-24-res.dat" index 0 using 1:7 title "Error real 2"