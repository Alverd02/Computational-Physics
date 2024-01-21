set ylabel "Energia [J]"
set xlabel "t [s]"
set xrange [0:7]
set key bmargin center horizontal box

set term png
set output "figE4.png"

plot  "apartat2a.dat" index 1 using 1:4 title"cas -","apartat2a.dat" index 1 using 1:7 title"cas +","apartat2a.dat" index 2 using 1:4 title"cas -","apartat2a.dat" index 2 using 1:7 title"cas +"