# tbcode

This repository contains the code files used for simulations concerning superconducting host systems with impurities embedded in or adsorbed on their surfaces. The code is split into the following parts:

### 📌 Superconducting host calculations

The BdG.f90 file, as well as the config.dat, hoppings.dat and basisvectors.dat configuration files comprise the main body of the program. The Bogolubov-de Gennes equations are first solved self-consistently for the non-superconducting system, in order to calculate its chemical potential, as well as the electron density on each lattice and atom site. Using the converged values, another self-consistency cycle is initiated, this time with the order parameter switched on, so that we can simulate a superconductor in 1, 2 or 3 dimensions. The converged values of this cycle correspond to the order parameter's value on each lattice and atom site.

### 📌 Band Structure

This is not part of the main program, however the Bands.f90 file, as well as the bandpoints.dat configuration file can be used in order to obtain band structure data along specific directions.

### 📌 Green functions and Impurity Problem

The Green.f90 file, along with the greenconfig.dat and the impurities.dat or the impvals.txt configuration files, calculates the Green function corresponding to the superconducting system we have modelled thus far and performs several calculations using this Green function. After automatically preparing a new configuration file, the BdGimp.f90 file is compiled in order to self-consistently solve Dyson's equation and calculate the Green function corresponding to the full system: the one containing the host superconductor, as well as the impurities.

### 📌 Miscellaneous

Other files that can be seen in this folder correspond to other procedures, such as the creation of specific impurity arrangements, or the calculation of the topological invariant (Pfaffian) for one-dimensional topological superconductors.
