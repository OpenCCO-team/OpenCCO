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
set output 'statTimeComp.pdf';
plot [][]  "statTimeVascu.dat" using 1:2 title "VascuSynth" w l lw 1, "statTimeCCO.dat" using 1:2 title "CCO 3D implicit domain" w l  lw 1 

