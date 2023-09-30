# script sencillo para plotear una funcion

# muestra los ejes
set xzeroaxis
set yzeroaxis

f(x)=sin(x)*exp(-abs(x))

set xrange[-4:4]
set yrange[-0.5:0.5]

set xlabel "x"
set ylabel "f(x)"
plot f(x)
pause -1
