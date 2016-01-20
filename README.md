# elastic - linear elasticity solver
Elastic is a computationally efficient linear elastic FEM solver for 2d and 3d simulations.

## Installation
clone this repository and build elastic:

` git clone https://github.com/ICStoolbox/LinearElasticity.git `

go to LinearElasticity directory

` cd LinearElasticity `

then create build directory and create Makefile

`mkdir build
cd build
cmake ..
make`

if no errors are produced, install the binary and library

` make install ` 

you can test the installation by entering the demos directory

` cd ../demos 
elastic demo1`

this should produce a result like

```- ELASTIC, Release 5.0c, Jan.19, 2016
   (C) Copyright 2006- , ICS-SU

 - LOADING DATA
    carre2.elas: 3 parameters
    carre2.mesh: 81 vertices, 32 edges, 128 triangles
 - COMPLETED: 0.001s

 ** MODULE ELASTIC: 5.0c
    Matrix and right-hand side assembly
    Solving linear system: 3.009659E-07 in 32 iterations
 ** COMPLETED: 0.001s

 - WRITING DATA
    carre2.sol: 81 data vectors
 - COMPLETED: 0.000s

 ** Cumulative time: 0.003s.
```
(Note: the version or release numbers may vary)

## Usage

## Documentation

## Platforms
elastic has been succesfully installed on Mac OSX and most Linux environments. 
   
## Authors & developers
elastic has been initiated by Maya de Buhan and Pascal Frey. Current team includes Charles Dapogny, Chiara Nardoni and Loic Norgeot.

Contributors to this project are warmly welcomed. You can help us to improve the code by 
* pull requests: 
* feature requests
* bug reports

Contact: 

## License and copyright
elastic is given under the [terms of the GNU Lesser General Public License](https://raw.githubusercontent.com/MmgTools/mmg/master/LICENSE).

Copyright © Université Pierre et Marie Curie, 2006 - .