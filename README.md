# Project Overview
This repository contains all analyses and manuscript files for the study “Climate vulnerable fish communities are also most exposed to human pressures”.

Authors: Fernanda Silva (fcads@aqua.dtu.dk), Daniël van Denderen, Antoni Vivó-Pons, Federico Maioli, Marcel Montanyès, Martin Lindegren

An archived version of the repository is available on Zenodo DOI

# Repository structure
R/ Contains all R scripts used for the data analysis. This folder includes code for data cleaning, randomization, statistical analysis, and outputs.

Data/ Contains the necessary data to generate the outputs of the manuscript.

Output/ Contains all graphics and figures resulted by the analysis and reported in the text. This folder serves as the central location for outputs referenced in the manuscript.

# Data availability
All datasets used in this project are publicly available from their respective sources:

Fish community’s data and biomass estimation compiled by FISHGLOB complemented with bottom-trawl survey data from Greenland and Iceland through the EU project B-USEFUL (https://b-useful.eu/), as well survey data from New Foundland kindly contributed by Fisheries and Oceans Canada, Newfoundland and Labrador Region (see reference in the manuscript).

Fish environmental tolerances available through AquaMaps (https://aquamaps.org/) and can be accessed using the R package “aquamapsdata”

Bottom temperature, salinity and oxygen data from the Copernicus Marine Environment Monitoring Service (CMEMS)

Human pressures from Halpern et al. 2019 (https://doi.org/10.1038/s41598-019-47201-9)

Conservation measures available at Protected Planet database through https://www.protectedplanet.net/en

# Reproducibility
To reproduce the analysis:

Clone or download this repository.
Open the .Rproj file in RStudio.
Install the required packages (if not already installed).

# Note
The processing and handling data is computationally intensive and may require substantial runtime. If your primary interest is the calculation of safety margins and modelling, you may start directly from code 5 - Proportion of safe sp.

Some datasets need to be downloaded from their respective website/repositories
