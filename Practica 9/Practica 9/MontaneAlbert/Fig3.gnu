set term png
set output "P9-23-24-fig3.png"

set title "T = 1280 ºC"

set key bmargin center horizontal box


set ylabel "T [ºC]"
set xlabel "iteracions"

plot  "res.dat" index 4 using 1:2 title"Jacobi" with lines  ,"res.dat" index 5 using 1:2 title"Sobrerelax w=1.45" with lines