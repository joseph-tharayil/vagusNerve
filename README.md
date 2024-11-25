# Simulation Insights on the Compound Action Potential

This repository provides the code necessary to replicate the results of the paper Simulation Insights on the Compound Action Potential. 

In broad terms, the workflow of the paper is as follows: 

First, finite element models of the vagus nerve are used, in Sim4Life, to calculate the E field produced by a stimulating electrode (two configurations are used in the paper), and the potential field produced by a virtual current applied between the two recording electrodes (several configurations are used in the paper). 

The stimulating E field is used to perform a neural titration, also in Sim4Life, for a sample of fibers of identical diameter. The results of this titration are used in the semi-analytic model to predict activation of the entire population of fibers in the vagus nerve. 

The potential field produced by the recording electrodes is interpolated along the center of each fascicle, and used in the semi-analytic model to calculate the eCAP, according to the reciprocity theorem. 

# Instructions

## System requirements

Our documentation and examples assume that you are running the semi-analytic models on a Linux system with slurm. The code to generate finite element models assumes that you are using the Sim4Life finite element platform (which is only compatible with Windows). This code has not been tested on any other system.

## Dependencies

The code for the semi-analytic models depends on MPI. The bash scripts provided in the *simulations* repository assumes that a module called *hpe-mpi* is available in an archive called *unstable*. This code also depends on several other python packages, which are automatically installed with setuptools when the package is installed.

## Installation
To install the code to run the semi-analytic models, create a virtual environment by running `python -m venv vagusEnv`. Then, simply run `pip install . ` to install this code in a virtual environment.

To run the finite element models, see the relevant documentation.

## Testing
The code to run the semi-analytic models can be tested by running `pytest tests`

## Replicating results from the paper

### Create finite element models
The finite element models used in this study can be created by extruding a 2D nerve histology cross-section (to which electrode geometries were manually added) along the longitudinal axis. Cross-sections are available on oSPARC for the model with the large recording electrode, with the small recording electrode on either the left or the right side, and for the horizontal and vertical stimulus electrode. They are in the form of iSeg project files; each cross-section consists of a .h5 file, a .prj file, and a .xmf file.

To create the nerve model, import these slices into the Sim4Life project files titles ExtrusionLines_*.smash (where * refers to each of the different models) and run the appropriate script in the create_fem_simulation/extrudeMesh folder to create the 2.5D model. Then run the scripts create_fem_simulation/createPerineuria/MakePatches.py and create_fem_simulation/createPerineuria/getFascicleDiameters.py (in that order) to create the thin layers for the electrode contacts and the perineuria. Note that these scripts will only produce the 2.5D mdoels used in teh EM simulations. For the neural simulations, you must use the provided FEM file, in which a region of the fascicles corresponding to the perinueria is segmented out. This ensures that the script which generates neurons does not place any within the perineuria.

###  Run finite element models
Once the finite element models are created in the previous step, they must be set up to run EM simulations to calculate the stimulating E field or the recording exposure function. These simulations are also available on oSPARC, with separate Sim4Life project files for the two stimulus cases (vertical and horizontal electrodes), for recording with the small electrodes on the left and on the right, and for bipolar recordings with the full cuff 1 cm from the stimulus electrode. Monopolar recordings and bipolar recordings 6cm from the stimulus electrode are all in the same Sim4Life project file (all_other_recording.smash). Simply run the simulations and export the data files in the postprocessing tab.

Alternatively, if you are setting up the simulations from scratch, perform the following steps to export the data:

For the EM simulations for recording exposure, run the scripts create_fem_simulation/extract_fem_results/create_splines_to_interpolate_over.py, then create_fem_simulation/extract_fem_results/InterpolatePhiFasc.py. These will save the exposure function to xlsx files, which are used by the semi-analytic models

For the EM simulations for stimulation, export the potential fields as a cache file. These will be used by the titration simulations.

### Run neuronal titrations
Sim4Life project files containing the neural titration are also available on oSPARC (titration_vertical_Rat.smash, titration_horizontal_Rat.smash, titration_vertical_Sundt.smash, titration_horizontal_Sundt.smash). Simply run the simulations and export the Titration Excel files in the postprocessing tab.

Alternatively, if you wish to set up the simulation yourself, run the script create_fem_simulation/createNeurons/FunctionalizeNeurons_Titrate_PerineuriaExcluded.py. Then import the cache file from the previous step and run the titration.

### Run semi-analytic models
Run all of the bash files in the *simulation* folder. Each bash file launches a python scripy that corresponds to a semi-analytic simulation of the eCAP, under various stimulus and recording parameters. In these python scripts, `stimulation_directory` refers to the paths to the Excel files output from the titration step. `recording_directory` refers to the paths to the Excel files containing the exposure functions. Note that you will have to change the paths in the bash files and corresponding python scripts to match the location of the input and output folders on your system.

### Create plots
Running all of the cells in the notebook provided will generate the majority of the figures in the paper. Note that you will need to change the hard-coded paths to match the location of the output folders on your system.
The other figure panels are created in Sim4Life. To create these panels, open the Sim4Life finite element models for the stimulation EM simulations for the two different electrode orientiations, and run the corresponding scripts in the figure_panels_s4l folder. These will change the colors of the fascicles.

# Citation
If you use this software, we kindly ask you to cite the following publication:
[Tharayil et al. Simulation Insights on the Compound Action Potential. *bioRxiv, (2024)*](https://www.biorxiv.org/content/10.1101/2024.10.16.618681v1.abstract)

# Acknowledgment
The development of this software was supported by funding to the Blue Brain Project, a research center of the École polytechnique fédérale de Lausanne (EPFL), from the Swiss government's ETH Board of the Swiss Federal Institutes of Technology.
 
Copyright (c) 2024 Blue Brain Project/EPFL
