set term png
set output 'P3-23-24-fig1.png'
set title 'Distancia al origen i su derivada en funcion de E'
set ylabel 'D(E) [U.A.]'
set xlabel 'E
plot 'P3-23-24-res.dat' index 1 using 1:2 title"D(E)",'P3-23-24-res.dat' index 1 using 1:3 title"D(E)'"