set term png
set output "P10-1819-fig2.png"

set title "Estacionari"

set key bmargin center horizontal box


set ylabel "T [ºC]"
set xlabel "t [s]"

plot  "apartat2a.dat"using 1:2 with lines title"x=6 cm","apartat2a.dat" using 1:3 with lines title"x=42 cm","apartat2a.dat" using 1:4 with lines title"x=126 cm"