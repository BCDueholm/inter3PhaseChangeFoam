This repository contains the solver and a benchmark case for OpenFOAM v9 from the Foundation and v2106 from ESI-OpenCFD.

The code and cases has been tested on debian 10, and issues might arise when trying to install it on Linux subsystems for windows.

The solver is an incompressible VOF solver capable of handling tree fluids with phase change between 1 and 2. The solver was developed to simulate cavitation in injectors.

## Prerequisite

To use the code and run the cases it is required to have either a working installation of.
- OpenFOAM v9
- OpenFOAM v2106
The postprocessing done in the banchmark cases can be done either using 
- python 3.8 with up to date versions of
	- numpy
	- matplotlib
- gnuplot 5.2 patchlevel 6

## Installation

To install the solver have one of the two versions of OpenFOAM sourced in a terminal and execute Allwmake in inter3PhaseChangeFoam in the respective version.

Further information can be found in the readme in the solver.

## Running the case

There is an Allrun script located in the benchmark case called phaseChangeColumn in the respective version, that will execute both part of the benchmark case.
Further information can be found in the readme in the benchmark case.
