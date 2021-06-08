# LiverCCO



Implementation of the CCO algorithm based on [[Schreiner etal. 93](https://github.com/kerautret/LiverCCO/blob/main/Refs/schreiner90.pdf)] algorithm and note of [[Clara Jaquet and Hugues Talbot  Feb 2021](https://github.com/kerautret/LiverCCO/blob/main/Refs/ccoJacquetHugues.pdf)]




## Kamiya’s algorithm for  local optimisation

The choosen implementation use the non linear solver of Ceres Solver: http://ceres-solver.org/installation.html#


Reproducting result of figure 8 (a) page 14 :
./testKamiya --x0 100 --y0 100 --x1 70 --y1 20 --x2 130 --y2 20 -s figure8_a.svg

Reproducting result of figure 8 (b) page 14 :
./testKamiya --x0 100 --y0 100 --x1 70 --y1 40 --x2 100 --y2 20 -s figure8_b.svg



Reproducting result of figure 8 (c) page 14 :
./testKamiya --x0 100 --y0 100 --x1 70 --y1 40 --x2 100 --y2 20 -s figure8_b.svg


Reproducting result of figure 8 (d) page 14 :
 ./testKamiya --x0 100 --y0 100 --x1 70 --y1 40 --x2 100 --y2 20 -s figure8_d.svg -R 0.25








   
   
