set ylabel "theta' [rad/s]"
set xlabel "t [s]"
set yrange [-0.1:0.15]

set term png
set output "P6-23-24-fig2.png"

f(x) = -0.02*sqrt(3.71/0.45)*sin(sqrt(3.71/0.45)*x)
plot  "P7-23-24-res.dat" index 0 using 1:2 title "Euler", f(x) title "theta' using sin(theta)=theta","P7-23-24-res.dat" index 0 using 1:3  title "Adams-Basford" 