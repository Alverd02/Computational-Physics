set ylabel "y(x)"
set xlabel "x"
set xzeroaxis
set yzeroaxis
set xrange[0:2]
set yrange[-2:5]


set term png
set output "P3-23-24-fig1.png"
plot (-57./160*pi + (57./80 + 17./20*pi)*x - (17./10 + pi/2.)*x**2 + x**3)*sinh(x) title"f(x)",((57./80 + 17./20*pi) - 2*(17./10 + pi/2.)*x + 3*x**2)*sin(x) + (-57./160*pi + (57./80 + 17./20*pi)*x - (17./10 + pi/2.)*x**2 + x**3)*cosh(x) title"f(x)'"