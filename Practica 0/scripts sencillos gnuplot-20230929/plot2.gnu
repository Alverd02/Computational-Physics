# script sencillo para plotear una funcion
# y comparala con un fichero de datos

# muestra los ejes
set xzeroaxis
set yzeroaxis

f(x)=sin(x)*exp(-abs(x))

set xrange[-4:4]
set yrange[-0.5:0.5]

set xlabel "x"
set ylabel "f(x)"
plot f(x), "datos.dat" u 1:2 t"Medidas experimenales" ps 2
pause -1
