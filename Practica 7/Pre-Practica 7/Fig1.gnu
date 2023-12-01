set ylabel "theta' [rad/s]"
set xlabel "t [s]"


set term png
set output "P6-23-24-fig1.png"


plot  "P7-23-24-res.dat" index 1 using 1:2 title "Euler","P7-23-24-res.dat" index 1 using 1:4 title "Adams-Bashford"