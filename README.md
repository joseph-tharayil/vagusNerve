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

### Create and run finite element models
Open the Sim4Life model containing the 2D cross-section of the nerve, and run the appropriate script in the create_fem_simulation/extrudeMesh folder to create the 2.5D model. Then run the scripts create_fem_simulation/createPerineuria/MakePatches.py and create_fem_simulation/createPerineuria/getFascicleDiameters.py (in that order) to create the thin layers for the electrode contacts and the perineuria. Note that these scripts will only produce the 2.5D mdoels used in teh EM simulations. For the neural simulations, you must use the provided FEM file, in which a region of the fascicles corresponding to the perinueria is segmented out. This ensures that the script which generates neurons does not place any within the perineuria.

For the EM simulations for recording exposure, run the scripts create_fem_simulation/extract_fem_results/create_splines_to_interpolate_over.py, then reate_fem_simulation/extract_fem_results/InterpolatePhiFasc.py. These will save the exposure function to xlsx files, which are used by the semi-analytic models

For the EM simulations for stimulation, export the potential fields as a cache file. These will be used by the titration simulations.

### Run neuronal titrations

For the titration simulations, run the script create_fem_simulation/createNeurons/FunctionalizeNeurons_Titrate_PerineuriaExcluded.py. Then import the cache file from the previous step and run the titration.

### Run semi-analytic models
Run all of the bash files in the *simulation* folder. Each bash file corresponds to a semi-analytic simulation of the eCAP, under various stimulus and recording parameters. 

### Create plots
Running all of the cells in the notebook provided will generate the majority of the figures in the paper.
The other figure panels are created in Sim4Life. To create these panels, open the Sim4Life finite element models for the stimulation EM simulations for the two different electrode orientiations, and run the corresponding scripts in the figure_panels_s4l folder. These will change the colors of the fascicles.

# Citation
If you use this software, we kindly ask you to cite the following publication:
Tharayil et al. Simulation Insights on the Compound Action Potential. *bioRxiv, (2024)*

# Acknowledgment
The development of this software was supported by funding to the Blue Brain Project, a research center of the École polytechnique fédérale de Lausanne (EPFL), from the Swiss government's ETH Board of the Swiss Federal Institutes of Technology.
 
Copyright (c) 2024 Blue Brain Project/EPFL
