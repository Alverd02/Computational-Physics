set term png
set output "P6-23-24-fig4.png"
set multiplot layout 2, 1  

set xlabel 't [s]'
set ylabel 'Energia [J]'
set key at graph 0.22,0.95
set yrange [-1000:1700]
plot 'P7-23-24-res.dat' index 2 using 1:2 title"Epot",'P7-23-24-res.dat' index 2 using 1:3 title"Etot"

set xlabel 't [s]'
set ylabel 'Energia [J]'
set yrange [-1000:10000]
plot 'P7-23-24-res.dat' index 2 using 1:4 title"Epot",'P7-23-24-res.dat' index 2 using 1:5 title"Etot"

unset multiplot