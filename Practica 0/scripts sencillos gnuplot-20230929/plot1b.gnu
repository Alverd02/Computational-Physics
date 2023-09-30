# script sencillo para plotear una funcion
# y generar un fichero png con la figura

set term png
set output "figuplot1b.png"

# muestra los ejes
set xzeroaxis
set yzeroaxis

f(x)=sin(x)*exp(-abs(x))

set xrange[-4:4]
set yrange[-0.5:0.5]

set xlabel "x"
set ylabel "f(x)"
plot f(x)
