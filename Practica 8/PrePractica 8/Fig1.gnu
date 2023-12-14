set term png
set output "P6-23-24-fig1.png"
set multiplot layout 2, 1  


set key bmargin center horizontal box


set ylabel "PHI"
set xlabel "x"

plot  "resultatsprepra8.dat" index 0 using 1:3 title"E1 400" with lines,"resultatsprepra8.dat" index 1 using 1:3 title"E2 400" with lines,"resultatsprepra8.dat" index 2 using 1:3 title"E3 400" with lines,"resultatsprepra8.dat" index 3 using 1:3 title"E4 400" with lines

set ylabel "PHI"
set xlabel "x"

plot  "resultatsprepra8.dat" index 4 using 1:3 title"E1 20" with lines,"resultatsprepra8.dat" index 5 using 1:3 title"E2 20" with lines,"resultatsprepra8.dat" index 6 using 1:3 title"E3 20" with lines,"resultatsprepra8.dat" index 7 using 1:3 title"E4 20" with lines

unset multiplot
