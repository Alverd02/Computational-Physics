set ylabel "theta' [rad/s]"
set xlabel "t [s]"
set key bmargin center horizontal box

set term png
set output "P6-23-24-fig2.png"


plot  "P7-23-24-res.dat" index 1 using 1:2 title "Euler","P7-23-24-res.dat" index 1 using 1:4 title "Euler Segon ordre"