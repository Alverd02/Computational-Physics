set term png
set output "P6-23-24-fig4.png"
set multiplot layout 2, 1  
set key bmargin center horizontal box

set xlabel 't [s]'
set ylabel 'Energia [J]'

plot 'P7-23-24-res.dat' index 2 using 1:2 title"Epot Euler" with lines,'P7-23-24-res.dat' index 2 using 1:3 title"Etot Euler" with lines,'P7-23-24-res.dat' index 2 using 1:6 title"Epot A-B" with lines,'P7-23-24-res.dat' index 2 using 1:7 title"Etot A-B" with lines

set xlabel 't [s]'
set ylabel 'Energia [J]'

plot 'P7-23-24-res.dat' index 2 using 1:4 title"Epot Euler" with lines,'P7-23-24-res.dat' index 2 using 1:5 title"Etot Euler" with lines,'P7-23-24-res.dat' index 2 using 1:8 title"Epot A-B" with lines,'P7-23-24-res.dat' index 2 using 1:9 title"Etot A-B" with lines

unset multiplot