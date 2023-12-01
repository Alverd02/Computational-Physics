set term png
set output "P6-23-24-fig4.png"
set key bmargin center horizontal box

set xlabel 't [s]'
set ylabel 'Energia [J]'

plot 'P7-23-24-res.dat' index 2 using 1:2 title"Etot Euler" with lines,'P7-23-24-res.dat' index 2 using 1:3 title"Ekin Euler" with lines,'P7-23-24-res.dat' index 2 using 1:4 title"Etot Euler 2" with lines,'P7-23-24-res.dat' index 2 using 1:5 title"Ekin Euler 2" with lines

