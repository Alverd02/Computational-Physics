set term png
set output "P10-1819-fig3.png"

set title "Evolució mitjana Ferro i Or"

set key bmargin center horizontal box


set ylabel "T [ºC]"
set xlabel "t [s]"

plot  "apartat2b.dat" index 0 using 1:2 with lines title"x=Ferro","apartat2b.dat" index 1 using 1:2 with lines title"Or"