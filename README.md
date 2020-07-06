# microbiome-hostspan

This repository includes the code used to produce the results presented in *Zapien, et al.*, 2020.

The aim is to study the effect of the host lifespan on a neutral model of the microbiome.

## Organization

1. Figures (**microbiome-hostspan/figures**).

Scripts to produce all the figures.

2. Numerics (**microbiome-hostspan/numerics**).

Scripts to solve the model numerically using the master equation.

3. Simulations (**microbiome-hostspan/simulations**).

Scripts to produce stochastic simulations.

## Usage

Scripts are structured in the following way:

`*_sc.py` contains the source code.

`*_exe.py` produces the output.

`*_par.py` contains user-defined parameters.

## Requirements

All the code has been developed in Python 3.6.
