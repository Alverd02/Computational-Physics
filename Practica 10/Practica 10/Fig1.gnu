set term png
set output "P10-1819-fig1.png"

set title "Estacionari"

set key bmargin center horizontal box


set ylabel "T [ºC]"
set xlabel "x [cm]"

plot  "apartat1.dat" index 0 using 1:2 with lines title"alpha = 0.00004","apartat1.dat" index 1 using 1:2 with lines title"alpha = 0.0003","apartat1.dat" index 2 using 1:2 with lines title"alpha = 0.0025"