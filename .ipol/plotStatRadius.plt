set terminal pdf
set output 'statRadiusBifLevel.pdf';


#set logscale
set grid;
#set grid xtics noxtics;
#set grid ytics noytics;

plot [][] "stat.dat" using 1:2:3 title "mean radius" w yerrorbars,  '' using 1:2 with lines  
 
#plot [] [2:1000] "forme15pt300max20NoMean.dat" using (20/exp($1)):(exp($2)) title "profile Mean" w lp,"forme15pt300max20NoMean.dat" using (20/exp($1)):(exp($4)) title "profile Min" w l, "forme15pt300max20NoMean.dat" using (20/exp($1)):(exp($5)) title "profile Max" w l,  x**(1.0/1.0)*4; 