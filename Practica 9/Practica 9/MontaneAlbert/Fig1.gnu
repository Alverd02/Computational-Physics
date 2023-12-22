set term png
set output "P9-23-24-fig1.png"

set title "T = 15 ºC"

set key bmargin center horizontal box


set ylabel "T [ºC]"
set xlabel "iteracions"

plot  "res.dat" index 0 using 1:2 title"Jacobi" with lines  ,"res.dat" index 1 using 1:2 title"Sobrerelax w=1.45" with lines