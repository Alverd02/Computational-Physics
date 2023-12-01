set xlabel "theta [rad]"
set ylabel "theta' [rad/s]"
set key bottom right

set term png
set output "P6-23-24-fig5.png"


plot  "P7-23-24-res.dat" index 3 using 3:2 title "2sqrt(h/l)+0.04","P7-23-24-res.dat" index 3 using 5:4 title "2sqrt(h/l)-0.04"