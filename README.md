# OpenCCO

[![CI (linux/macOS)](https://github.com/kerautret/LiverCCO/actions/workflows/build.yml/badge.svg)](https://github.com/kerautret/LiverCCO/actions/workflows/build.yml)



This is the Implementation of the CCO algorithm based on [[Schreiner etal. 93](https://github.com/kerautret/LiverCCO/blob/main/Refs/schreiner90.pdf)] algorithm and note of [[Clara Jaquet and Hugues Talbot  Feb 2021](https://github.com/kerautret/LiverCCO/blob/main/Refs/ccoJacquetHugues.pdf)]
and proposed to the IPOL journal: 
"OpenCCO: An Implementation of Constrained Constructive Optimization for Generating 2D and 3D Vascular Trees" 

version 1.0 IPOL article ( http://dx.doi.org/xxxx )

Copyright (c) 2023 by B. Kerautret;  Phuc Ngo, N. Passat H. Talbot and C. Jaquet
License: GNU LESSER GENERAL PUBLIC LICENSE (the "LGPL").
Version 1.0 04/04/2023


## Dependencies

- CMake >= 3.5
- [DGtal](https://github.com/DGtal-team/DGtal) 
   (eventually with QGLViewer to have 3D display tools but not mandatory)
- [DGtalTools-contrib](https://github.com/DGtal-team/DGtalTools-contrib.git) (optional) : provide **graphViewer** tool
- [DGtalTools](https://github.com/DGtal-team/DGtalTools) (optional) : provides **volBoundary2obj** and **meshViewer** tools that are used in the online demonstration.
- Ceres-seolver: can be installed from the source included in the ext directory or using the following:

Install dependencies CMake, ceres-solver Boost (for DGtal) :
```
sudo apt-get install cmake libboost-dev libceres-dev libceres1
```

## Installation


1. mkdir build && cd build
2. cmake .. 
3. make

Then the binary file "generateTree2D" and "generateTree3D" will be available in the "/build/bin" directory.


### run 
---------------
For more details see IPOL Journal article available here: 
 http://dx.doi.org/10.5201/ipol.xxx






# Changelog

Version 1.0 (03/04/2023):  initial online version


