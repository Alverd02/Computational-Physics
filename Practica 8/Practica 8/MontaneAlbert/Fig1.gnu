set term png
set output "P8-23-24-fig1.png"
 
set xrange [-7:1]

set key bmargin center horizontal box


set ylabel "PHI"
set xlabel "x"

plot  "integral.dat" index 0 using 1:2 title"E1" with lines,"integral.dat" index 1 using 1:2 title"E2" with lines,"integral.dat" index 2 using 1:2 title"E3" with lines,"integral.dat" index 3 using 1:2 title"E4" with lines

