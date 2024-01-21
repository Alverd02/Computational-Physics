set ylabel "Theta' [rad/s]"
set xlabel "Theta [rad]"

set key bmargin center horizontal box

set term png
set output "figE3.png"

plot  "apartat2a.dat" index 0 using 1:2 title"cas +","apartat2a.dat" index 0 using 3:4 title"cas -"