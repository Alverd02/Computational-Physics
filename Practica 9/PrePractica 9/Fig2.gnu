set term png
set output "P9-23-24-fig2.png"

set title "Sobrerelax w=1.5"

set key bmargin center horizontal box


set ylabel "T [ºC]"
set xlabel "x"

plot  "res.dat" index 1 using 1:2 title"T=2ºC" with lines,"res.dat" index 3 using 1:2  title"T=14ºC" with lines,"res.dat" index 5 using 1:2 title"T=230ºC" with lines 