# Project Name

Model code for: 

## Installation

Create conda environment using provided environment.yml file.

## Content

- Experimental-datafiles folder. Contains text files with experimental measurements.
  - growth_rates.txt: experimentally measured growth rate of all strains, in batch culture
  - : measured range expansion patterns (equilibrium frequency, sector size, overall growth)
- Processed-datafiles. Contains model output. All files are created by running code below.
  - community_data_mean.csv: processed format of community_data.csv (averaged over replicates)
  - fit_parameters.txt: fitted model parameters, output from fit_parameters_range_expansion.ipynb
  - model_predictions.csv: model predictions, output from plot_model_predictions_range_expansion.ipynb
- community.py: python code of community class, implements community predictions using analytical equations from [S1 Text from Ref 1.](https://doi.org/10.1371/journal.pcbi.1009877.s001)
- fit_parameters_range_expansion.ipynb: jupyter notebbok used to fit model parameters
- plot_model_predictions_range_expansion.ipynb: jupyter notebook used to make model predictions

## Usage

1. run fit_parameters_range_expansion.ipynb to refit model parameters
2. run plot_model_predictions_range_expansion.ipynb to recreate model predictions

[Ref 1]: van Vliet S, Hauert C, Fridberg K, Ackermann M, Dal Co A (2022) Global dynamics of microbial communities emerge from local interaction rules. PLOS Computational Biology 18(3): e1009877. https://doi.org/10.1371/journal.pcbi.1009877
