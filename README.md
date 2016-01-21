# elastic
Elastic is a finite element solver for linear elasticity problems in two and three dimensions.

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
After compiling elastic as described above, you should have an executable file in your $HOME/bin directory. If your PATH variable is correctly set to this directory, elastic can be called as follows:

    elastic mesh_file[.mesh]
 
A full description of all parameters that can be specified in the command line or in a parameter file [file.elas] can be found in the project [wiki](https://github.com/ICStoolbox/LinearElasticity/wiki).

#### Authors & contributors
* elastic has been initiated by Maya de Buhan (Université Paris Descartes) and Pascal Frey (Université Pierre et Marie Curie). Current team includes Charles Dapogny (Université Joseph Fourier), Chiara Nardoni (Université Pierre et Marie Curie) and Loic Norgeot (Université Pierre et Marie Curie).
* Contributors to this project are warmly welcomed. 

#### License
elastic is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).
