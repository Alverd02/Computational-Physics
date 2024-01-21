set ylabel "T [ºC]"
set xlabel "t [s] "


set term png
set output "figE3.png"


plot  "apartat2b.dat" index 0 using 2:1 title"Sin fuente","apartat2b.dat" index 1 using 2:1 title"Con fuente"