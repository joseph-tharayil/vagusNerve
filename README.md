# Simulation Insights on the Compound Action Potential

This repository provides the code necessary to run the semi-analytic model of the evoked compound action potential of the vagus nerve, as described in Simulation Insights on the Compound Action Potential. 

# Instructions

 ## System requirements

Our documentation and examples assume that you are running BlueRecording on a Linux system with slurm. This code has not been tested on any other system.

## Dependencies

This code depends on MPI. The bash scripts provided in the *simulations* repository assumes that a module called *hpe-mpi* is available in an archive called *unstable*. This code also depends on several other python packages, which are automatically installed with setuptools when the package is installed.

## Installation
Create a virtual environment by running `python -m venv vagusEnv`. Then, simply run `pip install . ` to install this code in a virtual environment.

## Testing
The code can be tested by running `pytest tests`

## Replicating results from the paper
Run all of the bash files in the *simulation* folder. Each bash file corresponds to a semi-analytic simulation of the eCAP, under various stimulus and recording parameters. Then run all of the cells in the notebook provided. This will generate the majority of the figures in the paper.
Some figure panels are created in Sim4Life.

# Citation
If you use this software, we kindly ask you to cite the following publication:
Tharayil et al. Simulation Insights on the Compound Action Potential. *bioRxiv, (2024)*

# Acknowledgment
The development of this software was supported by funding to the Blue Brain Project, a research center of the École polytechnique fédérale de Lausanne (EPFL), from the Swiss government's ETH Board of the Swiss Federal Institutes of Technology.
 
Copyright (c) 2024 Blue Brain Project/EPFL
