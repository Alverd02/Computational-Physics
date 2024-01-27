set term png
set output "Exa-jan-24-fig1.png"


set key bmargin center horizontal box


set ylabel "phi"
set xlabel "z"

plot  "Exa-jan-24-res1.dat" index 0 using 1:2 title"z0 = 0.3","Exa-jan-24-res1.dat" index 1 using 1:2 title"z0 = 0.9",