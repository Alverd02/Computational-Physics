set term png
set output "Exa-jan-24-fig2.png"


set key bmargin center horizontal box


set ylabel "I2"
set xlabel "N"

plot  "Exa-jan-24-res1.dat" index 5 using 1:2:3 with errorbars title""