# OpenCCO

[![CI (linux/macOS)](https://github.com/kerautret/LiverCCO/actions/workflows/build.yml/badge.svg)](https://github.com/kerautret/LiverCCO/actions/workflows/build.yml)



This is the Implementation of the CCO algorithm based on [Schreiner and Buxbaum 1993](https://ieeexplore.ieee.org/document/243413/) algorithm and note of [Jaquet and Huges 2021], the work is submitted to the IPOL journal: 
"OpenCCO: An Implementation of Constrained Constructive Optimization for Generating 2D and 3D Vascular Trees" 

The online demonstration is available [here](https://ipolcore.ipol.im/demo/clientApp/demo.html?id=5555531082026) 


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

Then the binary file `generateTree2D` and `generateTree3D` will be available in the `/build/bin` directory.

## Source code 

The code of OpenCCO is written in C++. It is composed of the following classes in the `src` directory:
* `CoronaryArteryTree` contains the structure of a vascular tree and functions to create a bifurcation
* `DomainController` contains the construction of the domain of a vascular tree
* `helpers/GeomHelpers` contains functions to generate randomly terminal points, check intersections, Kamyia optimisation, ...
* `helpers/ExpandTreeHelpers` contains functions to  generate a vascular tree
* `helpers/XMLHelpers` contains functions to export a vascular tree into XLM file.

## Program parameters 
### Options for both 2D and 3D cases
* number of terminal segments/ending points (-n | --nbTerm INT, default=1000)
* perfusion area / volume (-a | --aPerf FLOAT, default=20000)
* value of the gamma parameter (-g | -gamma FLOAT, default=3)
* minimal distance to border (-m | --minDistanceToBorder FLOAT, default=5) *Works only with option --organDomain
* organ domain using a mask image (-d | --organDomain TEXT)
* output the resulting geaph as xml file (-x | --exportXML TEXT)
* use a squared implicit domain instead a sphere (-s | --squaredDom, default=sphere) *Works only without option --organDomain

### Specific options in 2D
* initial position of root (-p | --posInit INT INT, default=image center)
* output the result into EPS format (-o | --outputEPS TEXT, default=result.eps) 
* output the result into SVG format (-e | --exportSVG TEXT, default=result.svg) 

### Specific options in 3D
* initial position of root (-p | --posInit INT INT INT, default=image center)
* output the 3D mesh into OFF file (-o | --outputName TEXT, default=result.off)
* output the 3D mesh into text file ( -e | --export TEXT)
* display 3D view using QGLViewer (--view)


## Typical 2D tree generation

### Generate vascalar tree on an implicit square domain 
(for example, 3000 ending points (-n) and a perfusion area of 20000 (-a) in a square domain (-s)):
```
./build/bin/generateTree2D -n 3000 -a 20000  -s
```
You will obtain such a display:

![Capture d’écran 2023-04-03 à 03 06 25](https://user-images.githubusercontent.com/772865/229390129-635e2863-5679-4065-b6c2-cefa921f79aa.png)


### Generate vascalar tree on the domain defined from the boudary of a binary shape:
(for example, 3000 ending points (-n) and a perfusion area of 20000 (-a) in pre-defined domain (-d) from the image Samples/shape3.pgm):
```
./build/bin/generateTree2D -n 3000 -a 20000  -d Samples/shape3.pgm 
```
<img width="628" alt="Capture d’écran 2023-04-03 à 03 18 19" src="https://user-images.githubusercontent.com/772865/229391012-b9b6efc4-aa10-48de-ac5a-063d023f083f.png">



## Typical 3D tree generation

### Generate vascalar tree on the domain defined from the boudary of the bunny.obj:
The following commands permits to generate a vascular tree starting with a specific 3D point: 
(for example, 3000 ending points (-n) and a perfusion volume of 20000 (-a) in pre-defined domain (-d) from the file Samples/bunnyThickBdr.vol and with the initial position of root (-p) at (143, -107, 7)):
```
 ./build/bin/generateTree3D -n 3000 -a 20000  -d Samples/bunnyThickBdr.vol --view -m 1 -p 143 -107 7
 ```
 
<img width="616" alt="Capture d’écran 2023-04-03 à 02 46 55" src="https://user-images.githubusercontent.com/772865/229388906-2035b721-f4f6-4f9c-bb2b-1490b5a86187.png">


### Generate vascalar tree on an implicit square domain :
(for example, 3000 ending points (-n) and a perfusion area of 20000 (-a) in a square domain (-s)):
```
./build/bin/generateTree3D -n 3000 -a 20000  -s meshViewer result.off
```
(or used directly the --view option)
You will obtain such type of visualisation:
<img width="912" alt="Capture d’écran 2023-04-03 à 02 58 25" src="https://user-images.githubusercontent.com/772865/229389571-ccac9ca2-a560-4b1b-acce-7ce9d825efa5.png">

For more details see IPOL Journal article available here: 
 http://dx.doi.org/10.5201/ipol.xxx

# Acknowledgements
This work was supported by the French Agence Nationale de la Recherche (grants ANR-18-CE45- 0018, ANR-18-CE45-0014, ANR-20-CE45-0011).

# Changelog
Version 1.0 (03/04/2023):  initial online version


