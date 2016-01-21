# elastic - linear elasticity solver
Elastic is a computationally efficient linear elastic FEM solver for 2d and 3d simulations.

#### Installation
1. you will need to install the [ICS Commons Library](https://github.com/ICStoolbox/Commons) on your system. 
Please refer to the instructions provided on the ICS Commons Library page in order to install this library.

2. clone this repository and build elastic:

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

you can test the installation and look at examples by entering the [demos](demos) directory.

#### Usage
elastic can be used in two ways:
* as a standalone binary code. usage: elastic file.mesh
* as a library that can be plugged into C/C++ codes.
See the project [wiki](https://github.com/ICStoolbox/LinearElasticity/wiki) for more details.

#### Authors & contributors
* elastic has been initiated by Maya de Buhan (Université Paris Descartes) and Pascal Frey (Université Pierre et Marie Curie). Current team includes Charles Dapogny (Université Joseph Fourier), Chiara Nardoni (Université Pierre et Marie Curie) and Loic Norgeot (Université Pierre et Marie Curie).
* Contributors to this project are warmly welcomed. 

#### License
elastic is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).

Copyright © Université Pierre et Marie Curie, 2006 - .
