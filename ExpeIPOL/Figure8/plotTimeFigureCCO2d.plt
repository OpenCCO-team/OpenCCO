set xlabel "Number of vertex (k) " 
set ylabel "User execution time (sec)"
set grid;

set size ratio 0.8
set bars small
set lmargin at screen 0.13;
set rmargin at screen -1.29;
set bmargin at screen 0.15;
set tmargin at screen 0.95;

set terminal pdf size 3.9,3;
set output 'statTimeCompCCO2d.pdf';
plot [][]  "statTimeCCO2DimplicitDom.dat" using 1:2 title "CCO 2D with masked domain" w l lw 1, "statTimeCCO2Dimplicit.dat" using 1:2 title "CCO 2D implicit domain" w l  lw 1 

