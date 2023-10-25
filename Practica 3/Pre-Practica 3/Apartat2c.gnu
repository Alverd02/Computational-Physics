set term png
set output 'P3-23-24-fig2.png'
set title 'Convergencia del método de Newton-Raphson'
set ylabel 'Valor de la raiz'
set xlabel 'Número de Iteraciones'
plot 'P3-23-24-res.dat' index 1 using 1:2 title"E0=0.2",'P3-23-24-res.dat' index 2 using 1:2 title"E0=0.65",'P3-23-24-res.dat' index 5 using 1:2 title"E0=2.4"