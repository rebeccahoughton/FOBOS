# FOBOS: Few Observation Binary Orbit Solver

Note: This project is still being edited.

## About the project

FOBOS is an open source software written in Fortran 90, which takes 2-4 astrometric observations of a binary or triple star/exoplanetary system and a text tile which can be used to produce the probability distribution function of the system's orbital parameters. 

There are several versions of FOBOS available:
* FOBOS-binaries.f90 - Designed for use on a binary system with two epochs of astrometric data
* FOBOS-triples.f90 - Designed for use on a triple system with two epochs of astrometric data
* FOBOS-ME4.f90 - Designed for use on a binary system with four epochs of astrometric data

## Getting started

The instructions below detail how to set up this project to run locally. 

To make any changes to the source code files (for example, the number of threads for the simulation to run on, the output file directory, etc.), the code will need to be recompiled. FOBOS was designed for use with the gfortran compiler with OpenMP support. Alternative fortran compilers can be used but the following instructions may not work. 

Compile the code with the command:
~~~
  gfortran -fopenmp FOBOS-binaries.f90 -o run.exe
~~~
Where run.exe is the executable used to run the simulation. 

## Output

The simulation will output a file called params-'id' where 'id' is a user defined identifier for the system. 

## License

Distributed under the MIT License. See LICENSE.txt for more information

## Contact

Rebecca Houghton - @astro_rebecca - rhoughton1@sheffield.ac.uk
Project link: https://github.com/rebeccahoughton/FOBOS
