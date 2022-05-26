# FOBOS: Few Observation Binary Orbit Solver

Note: This project is still being edited.

## About the project

FOBOS is an open source orbital fitting software written in Fortran 90. It takes 2-4 epochs of relative astrometry for a binary or triple star/exoplanetary system and produces a text file of data which can be used to generate probability distribution functions for the system's orbital parameters. 

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

To run the simulation, you will also need the following files in the same directory:
* input.txt - contains the epoch and relative astronomy data in the same format as shown in the 'example' folder
* mass-dist.txt - contains estimates of the masses and distance of the objects/system
* sys-name.txt - contains the name of the system/identifier which will be used to name the output files.

## Output

The simulation will output a file called params-'id' where 'id' is a user defined identifier for the system. This file will contain 9 lines of metadata about the simulation. Following this, when using FOBOS-binaries.f90 and FOBOS-ME4.f90, there will be 5 columns with each row corresponding to the semi-major axis, eccentricity, inclination, orientation, and mean anomoly (left to right) of a synthetic system which fits the astrometric data.
When using FOBOS-triples.f90, the output file will contain 10 columns, starting with the semi-major axis of the secondary, semi-major axis of the tertiary, eccentricity of the secondary, etc. (left to right). 

## License

Distributed under the MIT License. See LICENSE.txt for more information

## Note

I am currently working on wrapping the FOBOS code in python to make it more accessible. I am more than happy to accept data and run it through FOBOS myself, and return the results as probability distributions rather than raw text files. Please see my contact details below. 

## Contact

Rebecca Houghton - @astro_rebecca - rhoughton1@sheffield.ac.uk
Project link: https://github.com/rebeccahoughton/FOBOS
