# TFG
#### Miquel Limón Vallés

This repository contains the code written for the final project of my degree in Physics.

### Main goals
The main goal of this project is to obtain a μ/U phase diagram of the 1D Fermi-Hubbard model using tensor network formalism. Using Julia and the libraries ITensors and ITensorMPS, I have made a set of functions to solve an iMPS-iTEBD problem using ITensor objects. Once this is polished and simulations can be run comfortably, i expect to obtain physically relevant results.

### Repository basics
The relevant files included here (ignoring test files) are the following:

- "*functions.jl*": contains a detailed collection of functions, the need for which arises as I get deeper into the project. These functions serve as an interface to work with ITensor-based iMPS comfortably, despite the official ITensor library not offering formal support for iMPS yet. 
- "*manual_itebd.jl*": using the functions mentioned above, I am trying to implement a simple infinite Time Evolving Block Decimation algorithm. My current project is simulating a single Up fermion in a 1D lattice with periodic boundary conditions and seeing that the expected solution is accomplished.
