set term png
set output "P9-23-24-fig2.png"

set title "T = 220 ºC"

set key bmargin center horizontal box


set ylabel "T [ºC]"
set xlabel "iteracions"

plot  "res.dat" index 2 using 1:2 title"Jacobi" with lines  ,"res.dat" index 3 using 1:2 title"Sobrerelax w=1.45" with lines