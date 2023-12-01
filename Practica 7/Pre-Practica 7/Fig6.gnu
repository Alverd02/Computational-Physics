set term png
set output "P6-23-24-fig6.png"


set xlabel 't [s]'
set ylabel 'Energia [J]'
set key left top
plot 'P7-23-24-res.dat' index 4 using 1:2 with lines title "Etot n=300",'P7-23-24-res.dat' index 5 using 1:2 with lines title "Etot n=1000",'P7-23-24-res.dat' index 6 using 1:2 with lines title "Etot n=2200",'P7-23-24-res.dat' index 7 using 1:2 with lines title "Etot n=14500"