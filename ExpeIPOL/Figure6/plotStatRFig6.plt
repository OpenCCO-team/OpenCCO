set xlabel "bifurcation level " 
set ylabel "mean segment diameter (mm)"
set grid;

set size ratio 0.8
set bars small
set lmargin at screen 0.13;
set rmargin at screen -1.29;
set bmargin at screen 0.15;
set tmargin at screen 0.95;

set terminal pdf size 3.9,3;
set output 'statRFig6.pdf';
plot [][]  "vascuStatR.dat" using 1:($2*25):($3*25) title "VascuSynth" w yerrorlines pt 0 lw 1, "ccoImplStatR.dat" using 1:2:3 title "CCO implementation" w yerrorlines pt 0 lw 1 

