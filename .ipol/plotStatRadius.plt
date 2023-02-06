set xlabel "bifurcation level \n (i.e number of proximal bifurcations)"
set ylabel "mean segment diameter (mm)"
set grid;


set terminal png
set output 'statRadiusBifLevel.png';
plot [][] "stat.dat" using 1:2:3 title "CCO algorithm" w yerrorlines


set terminal pdf
set output 'statRadiusBifLevel.pdf';
plot [][] "stat.dat" using 1:2:3 title "CCO algorithm" w yerrorlines

