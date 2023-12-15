set term png
set output "P8-23-24-fig2.png"
 


set key bmargin center horizontal box


set ylabel "E"
set xlabel "i"

plot  "convergencia.dat" index 0 using 1:2 title"E1" with lines,"convergencia.dat" index 1 using 1:2 title"E2" with lines,"convergencia.dat" index 2 using 1:2 title"E3" with lines,"convergencia.dat" index 3 using 1:2 title"E4" with lines

