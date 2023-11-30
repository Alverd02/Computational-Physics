set ylabel "theta [rad]"
set xlabel "theta' [rad/s]"


set term png
set output "P6-23-24-fig3.png"


plot  "P7-23-24-res.dat" index 1 using 3:2 title "Euler"