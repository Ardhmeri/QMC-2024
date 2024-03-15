Quantum Monte Carlo Simulation Repository
This repository contains Fortran code for performing Quantum Monte Carlo (QMC) simulations. QMC is a powerful computational method used to study quantum systems, particularly in the field of condensed matter physics and quantum chemistry.

Contents
Main program :
QuantumMonteCarlo.f90: Fortran source code for the general Quantum Monte Carlo algorithm.

Functions:
kinetic.f90: Fortran code for calculating the kinetic energy.
f_potential.f90: Fortran code for calculating the potential energy.
f_psi.f90: Fortran code for calculating the wave function.
local_energy.f90: Fortran code for computing the local energy.

Subroutines:
PuredifusionMonteCarlo.f90: Fortran  code for the Pure Diffusion Monte Carlo (PDMC) algorithm.
random_gauss.f90: Fortran code for generating Gaussian random numbers.
s_ave_error.f90: Fortran code for calculating the average error.
s_drift.f90: Fortran source code for implementing drift in the simulation.

Input data files for specific quantum systems:
H2+_ion.dat, H2_molecule.dat, H3+_ion.dat, H_atom.dat, He_atom.dat

Objectfiles:
 QuantumMonteCarlo.o , f kinetic.o ,f_potential.o , f psi.o , local_energy.o PuredifusionMonteCarlo.o , random_gauss.o , s_ave_error.o ,s_drift.o

qmc.x: Executable file generated after compiling the code.

To copile the code :
Use the Fortran compiler (gfortran, ifort, etc.) to compile each Fortran source code file (*.f90) into its corresponding object file (*.o).
                gfortran -c  QuantumMonteCarlo.f90 -o  QuantumMonteCarlo.o 
 Link all object files into the final executable
                  gfortran *.o -o qmc.x 
Run the compiled program:
                  ./qmc.x
 Afer you compile , the program will ask you to choose a number for a corresponding Input file , and depending on
 the number you choose , it will show you the result for that specific Quantum System .
 
                  
                  
                  
