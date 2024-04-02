# Overview

This repository contains the solver and a benchmark case for OpenFOAM v9 from the Foundation and v2106 from ESI-OpenCFD.

The code and cases have been tested on Debian 10.

The solver is an incompressible VOF solver capable of handling tree fluids with phase change between 1 and 2. The solver was developed to simulate cavitation in injectors.

# Prerequisites

To use the code and run the cases (including postprocessing), the following prerequisites are required:
- OpenFOAM 9 or OpenFOAM v2106
- Python 3.8 with up-to-date versions of
	- NumPy
	- Matplotlib
- Gnuplot 5.2 patchlevel 6

# Installation

To install the solver, execute `./Allwmake` in the respective `inter3PhaseChangeFoam` directory. For example, when using OpenFOAM 9, execute `cd v9/inter3PhaseChangeFoam && ./Allwmake`.

Further information can be found in the `readme` of the respective solver directories.

# Running the case

The cases can be run using the included `Allrun` scripts. For example, when using OpenFOAM 9, execute `cd v9/phaseChangeColumn && Allrun`.

Further information can be found in the `readme` of the respective benchmark cases.
