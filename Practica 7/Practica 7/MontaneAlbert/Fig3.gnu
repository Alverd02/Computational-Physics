set ylabel "theta' [rad/s]"
set xlabel "theta [rad]"
set key top left

set term png
set output "P6-23-24-fig3.png"


plot  "P7-23-24-res.dat" index 1 using 2:3 title "Euler","P7-23-24-res.dat" index 1 using 4:5 title "Euler Segon ordre"