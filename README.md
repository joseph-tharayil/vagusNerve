# Simulation Insights on the Compound Action Potential

This repository provides the code necessary to replicate the results of the paper [Simulation Insights on the Compound Action Potential](https://www.biorxiv.org/content/10.1101/2024.10.16.618681v1.full). 

In broad terms, the workflow of the paper is as follows: 

First, finite element models of the vagus nerve are used, in Sim4Life, to calculate the E field produced by a stimulating electrode (two configurations are used in the paper), and the potential field produced by a virtual current applied between the two recording electrodes (several configurations are used in the paper). 

The stimulating E field is used to perform a neural titration, also in Sim4Life, for a sample of fibers of identical diameter. The results of this titration are used in the semi-analytic model to predict activation of the entire population of fibers in the vagus nerve. 

The potential field produced by the recording electrodes is interpolated along the center of each fascicle, and used in the semi-analytic model to calculate the eCAP, according to the reciprocity theorem. 

# Instructions

## System requirements
The code to generate semi-analytic models has been confirmed to work on both Linux and Windows, with Python 3.13. It should also work on MacOS, and with other versions of python. The code to generate finite element models assumes that you are using the Sim4Life finite element platform (which is only compatible with Windows). This code has not been tested on any other system.

## Installation
To install the code to run the semi-analytic models, create a virtual environment by running `python -m venv vagusEnv`. Then, simply run `pip install . ` to install this code in a virtual environment.

This code requires finite element models which are available on [Zenodo](https://zenodo.org/records/15848037) or [oSparc](https://osparc.io/#/study/e03935d2-79a5-11ef-8269-0242ac174ae0).

To run the finite element models, see the relevant documentation.

## Testing
The code to run the semi-analytic models can be tested by running `pytest tests`

## Replicating results from the paper

### Create finite element models
The finite element models used in this study can be created by extruding a 2D nerve histology cross-section (to which electrode geometries were manually added) along the longitudinal axis. Cross-sections are available on oSPARC or Zenodo for the model with the large recording electrode, with the small recording electrode on either the left or the right side, and for the horizontal and vertical stimulus electrode. They are in the form of iSeg project files; each cross-section consists of a .h5 file, a .prj file, and a .xmf file.

To create the nerve model, import these slices into the Sim4Life project files titles `ExtrusionLines_*.smash` (where * refers to each of the different models) and run the appropriate script in the `create_fem_simulation/extrudeMesh` folder to create the 2.5D model. Then run the scripts `create_fem_simulation/createPerineuria/MakePatches.py` and `create_fem_simulation/createPerineuria/getFascicleDiameters.py` (in that order) to create the thin layers for the electrode contacts and the perineuria. Note that these scripts will only produce the 2.5D mdoels used in teh EM simulations. For the neural simulations, you must use the provided FEM file, in which a region of the fascicles corresponding to the perinueria is segmented out. This ensures that the script which generates neurons does not place any within the perineuria.

###  Run finite element models
Once the finite element models are created in the previous step, they must be set up to run EM simulations to calculate the stimulating E field or the recording exposure function. These simulations are also available on oSPARC and Zenodo, with separate Sim4Life project files (`Stimulus_Horizontal_v2.smash` and `Stimulus_Vertical.smash`) for the two stimulus cases (vertical and horizontal electrodes), and for recording with the small electrodes on the left and on the right (`SmallRecordingElectrode.smash` and `SmallRecordingElectorde_Otherside.smash`). Monopolar recordings and bipolar recordings 6cm from the stimulus electrode are all in the same Sim4Life project file (`all_other_recording.smash`). Simply run the simulations and export the data files in the postprocessing tab.

Alternatively, if you are setting up the simulations from scratch, perform the following steps to export the data:

For the EM simulations for recording exposure, run the `scripts create_fem_simulation/extract_fem_results/create_splines_to_interpolate_over.py`, then `create_fem_simulation/extract_fem_results/InterpolatePhiFasc.py`. These will save the exposure function to xlsx files, which are used by the semi-analytic models

For the EM simulations for stimulation, export the potential fields as a cache file. These will be used by the titration simulations.

### Run neuronal titrations
Sim4Life project files containing the neural titration are also available on oSPARC (`titration_vertical_Rat.smash`, `titration_horizontal_Rat.smash`, `titration_vertical_Sundt.smash`, `titration_horizontal_Sundt.smash`). Simply run the simulations and export the Titration Excel files in the postprocessing tab.

Alternatively, if you wish to set up the simulation yourself, run the script `create_fem_simulation/createNeurons/FunctionalizeNeurons_Titrate_PerineuriaExcluded.py`. Then import the cache file from the previous step and run the titration.

### Run semi-analytic models
Run all of the python scripts files in the *simulation* folder. Each corresponds to a semi-analytic simulation of the eCAP, under various stimulus and recording parameters. In these python scripts, `stimulation_directory` refers to the paths to the Excel files output from the titration step. `recording_directory` refers to the paths to the Excel files containing the exposure functions. Note that you will have to change the paths in the python scripts to match the location of the input and output folders on your system.

### Create plots
Running all of the cells in the notebook provided will generate the majority of the figures in the paper. Note that you will need to change the hard-coded paths to match the location of the output folders on your system.
The other figure panels are created in Sim4Life. To create these panels, open the Sim4Life finite element models for the stimulation EM simulations for the two different electrode orientiations, and run the corresponding scripts in the figure_panels_s4l folder. These will change the colors of the fascicles.

To generate the plots in Figure 3 in our paper, see [here](https://github.com/joseph-tharayil/vagusNerve/tree/main/validation). To generate the plots in Figure 4, see [here](https://github.com/joseph-tharayil/vagusNerve/blob/main/notebooks/ReadMe.md). To generate the plots in Figure 12, see [this repo](github.com/joseph-tharayil/vagusOptimization).

# Citation
If you use this software, we kindly ask you to cite the following publication:
[Tharayil et al. Simulation Insights on the Compound Action Potential. *PLOS Computational Biology, (2025)*](https://doi.org/10.1371/journal.pcbi.1013452)

# Acknowledgment
The development of this software was supported by funding to the Blue Brain Project, a research center of the École polytechnique fédérale de Lausanne (EPFL), from the Swiss government's ETH Board of the Swiss Federal Institutes of Technology.
 
Copyright (c) 2024 Blue Brain Project/EPFL
