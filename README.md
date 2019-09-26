## PGOPT

Parallel global optimization of gas phase and surface systems.

Copyright (C) 2018 Huanchen Zhai

Zhai, Huanchen, and Anastassia N. Alexandrova. "Ensemble-average representation of Pt clusters in conditions of catalysis accessed through GPU accelerated deep neural network fitting global optimization." *Journal of chemical theory and computation* **12** (2016): 6213-6226.

Zhai, Huanchen, and Anastassia N. Alexandrova. "Local Fluxionality of Surface-Deposited Cluster Catalysts: The Case of Pt7 on Al2O3." *The journal of physical chemistry letters* **9** (2018): 1696-1702.

## Basic User Guide

# Installation

This code contains three sub-packages. `ACNN` contains all core algorithms. `PGOPT` includes settings for supercomputer environment. `STMOLE` defines the interface to `VASP` and `TURBOMOLE`.

If you want to test the code in a non-supercomputer environment, `PGOPT` and `STMOLE` are not needed.

Here we will explain installation steps using a `docker` (https://docs.docker.com/) image of `anaconda2` (https://hub.docker.com/r/continuumio/anaconda). Please first install the `docker` application in your system, then open a terminal to continue. (You can also use `anaconda` with `python2` installed in your local environment without `docker`.) The following commands will download a image of `anaconda2`:

```
docker pull continuumio/anaconda
docker run -it --rm continuumio/anaconda /bin/bash
```

Now you are inside the `docker` container. We need to add the following `python` packages and several packages needed for compiling `fortran` code:

```
pip install theano reportlab dill
apt-get update
apt-get -y install gfortran g++ make vim
```

Now get the copy of `PGOPT`:

```
cd ~
git clone https://github.com/hczhai/PGOPT.git PGOPT-PROGRAMS
```

Compile the `fortran` code. This will generate some warnings but the process should not produce any error.

```
cd ~/PGOPT-PROGRAMS/ACNN/formod
make
```

Now add the following to `~/.bashrc` (you may need `vim ~/.bashrc` first):

```
BASE=~/PGOPT-PROGRAMS
export STMOLE_HOME=$BASE/STMOLE
export ACNNHOME=$BASE/ACNN
export PGOPTHOME=$BASE/PGOPT
export PATH=$STMOLE_HOME:$PATH
export PATH=$PGOPTHOME:$ACNNHOME:$PATH
```

Then apply the these environment changes:
```
source ~/.bashrc
```

# Gas Phase Cluster Generation

The main program is called `acnnmain` which can be found under `$ACNNHOME` defined previously. You need to prepare an input file for generating structures. There are some example input files under `ACNN/tests`. Here as an example, we will try to generate some gas phase Pt<sub>7</sub> structures using S-BLDA.

```
cd $ACNNHOME/tests/structure-generation
acnnmain pt7-gas.json
```

Note that the output directory is indicated in the input file. Here we can find the results in `./OUT-pt7-gas/fil_structs.xyz.0`. The structures are written in `XYZ` format. A single file will contains more than one structures. Visualization software such as `jmol` is useful for examine all structures within only one file. For example, if `jmol` is installed, you can type:

```
jmol ./OUT-pt7-gas/fil_structs.xyz.0
```

to look at the structures.

The `PGOPT` code also has its own visualization implementation. It will generate a PDF file containing images of all structures in the given input file. To generate the PDF, use the following input file (this only works after you run `acnnmain pt7-gas.json`):

```
acnnmain pt7-gas-draw.json
```

Then you will find the PDF in `./OUT-pt7-gas/report.pdf`.

# Surface Supported Cluster Generation

The following command will generate Pt<sub>7</sub> structures on alpha-Al<sub>2</sub>O<sub>3</sub> surface. The surface is described by a `XYZ` file. There are some example surface files under `$ACNNHOME/tests/surfaces`. The computational cell information is written in the comment line of the `XYZ` file. It can be either 3 numbers or 5 numbers. The program always assume the `Z` direction is normal to the surface plane. If the cell size is described by 3 numbers `n1 n2 n3`, then the `XYZ` components of cell axes are `a = (n1, 0, 0), b = (0, n2, 0), c = (0, 0, n3)`. If the cell size is described by 5 numbers `n1 n2 n3 n4 n5`, then the `XYZ` components of cell axes are `a = (n1, n2, 0), b = (n3, n4, 0), c = (0, 0, n5)`. Usually `n5` is larger than the actual height of the surface because of the added vacuum gap. So in the end of comment line there is an additional number in parenthesis, indicating the unit cell height along `Z` without vacuum gap.

The surface group and unit cell (the minimal unit cell, not the computational unit cell) information, indicated in the input file may help determining structure duplicates. If these information is unavailable, use the computational cell for unit cell and "P 1" for space group, and `[0.0, 0.0]` for space group transformation reference point.

```
acnnmain pt7-alpha.json
```

Then we can find the results in `./OUT-pt7-alpha/fil_structs.xyz.0`.

# Structure Filtering

The following command will try to find unique structures from a given example `XYZ` file containing some local minima (`$ACNNHOME/tests/data/pt4b4-local.xyz`).

```
cd $ACNNHOME/tests/filtering
acnnmain pt4b4-local.json
```

The unique structures will be in `./OUT-pt4b4-filter/fil_structs.xyz.0`. The additional file `./OUT-pt4b4-filter/fil_list.txt.0` shows how many duplicates of each unique structure appear in the original input `XYZ` file (the `multi` column). The other additional file `./OUT-pt4b4-filter/fil_corr.txt.0` lists the structural difference data. In this file, if one line ends with `*`, then the structure is selected as unique structure, because its structural difference to all previous unique structures is higher than the threshold. The second last column `mindm` shows the minimal value over structural difference to all previous structures.

If the structure filtering should be performed on surface support clusters, the `creation-surface` section should be given in input file, which contains the same information as that in the input file for creation.

# Neural Network Fitting

Note that Neural Network Fitting is only implemented for gas phase clusters containing only one type of element. Other research groups has published more general codes on this topic. For example, see `JCTC, 14(7), 2018, 3933-3942`.

The following command will try to fit a neural network based on an example Pt<sub>9</sub> data (`$ACNNHOME/tests/data/pt9-structs.xyz`). For realistic results, you need to change `sample_number` parameter from `[5000, 500, 500]` to `[200000, 20000, 20000]`. But then it will take a very long time (you may need a gpu). If there is no gpu to use, it will automatically switch to cpu.

```
cd $ACNNHOME/tests/nn_fitting
acnnmain pt9-fit.json
```

The fitted network will be in `./OUT-pt9-nn/fit_network.dill.0`.

Next we need to create some new structures, then use the network to optimize them.

```
acnnmain pt9-create.json
acnnmain pt9-opt.json
```

The optimized structures and energies will be in `./OUT-pt9-nn/opt_structs.xyz.0`.
