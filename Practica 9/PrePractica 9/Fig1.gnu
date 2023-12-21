set term png
set output "P9-23-24-fig1.png"

set title "Jacobi"

set key bmargin center horizontal box


set ylabel "T [ºC]"
set xlabel "iteracions"

plot  "res.dat" index 0 using 1:2 title"T=2ºC" with lines  ,"res.dat" index 2 using 1:2 title"T=14ºC" with lines,"res.dat" index 4 using 1:2 title"T=230ºC" with lines 
