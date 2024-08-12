# Auxotroph Range Expansion Model

Model code for: Engineering microbial consortia: uptake and leakage rate differentially shape community arrangement and composition.    
Estelle Pignon [1], Gábor Holló [1], Théodora Steiner [1], Simon van Vliet [1,2], Yolanda Schaerli [1]  
[1] Department of Fundamental Microbiology, University of Lausanne.  
[2] Biozentrum, University of Basel.

## Installation

Create conda environment using provided environment.yml file.
i.e. use `conda env create -f environment.yml` (install time several minutes)

## Content

- Experimental-datafiles folder. Contains text files with experimental measurements.
  - growth_rates.txt: experimentally measured growth rate of all strains, in batch culture
  - community_data.csv: measurements of range expansion patterns (equilibrium frequency, sector size, overall growth)
- Processed-datafiles folder. Contains model output. All files are created by running code below.
  - community_data_mean.csv: processed format of data in community_data.csv (data averaged over replicates), output from fit_parameters_range_expansion.ipynb
  - fit_parameters.txt: fitted model parameters, output from fit_parameters_range_expansion.ipynb
  - model_predictions.csv: model predictions, output from plot_model_predictions_range_expansion.ipynb
- community.py: python code of community class, implements community predictions using analytical equations from [S1 Text](https://doi.org/10.1371/journal.pcbi.1009877.s001) and [S2 text](https://doi.org/10.1371/journal.pcbi.1009877.s002) from [Ref 1](https://doi.org/10.1371/journal.pcbi.1009877).
- fit_parameters_range_expansion.ipynb: jupyter notebook used to fit model parameters
- plot_model_predictions_range_expansion.ipynb: jupyter notebook used to make model predictions

## Usage

1. run fit_parameters_range_expansion.ipynb to refit model parameters (runtime ~1min), saves fitted parameters to Processed-datafiles/fit_parameters.txt.
2. run plot_model_predictions_range_expansion.ipynb to recreate model predictions (runtime ~1min), saves figures to the Figures subfolder.

[Ref 1]: van Vliet S, Hauert C, Fridberg K, Ackermann M, Dal Co A (2022) Global dynamics of microbial communities emerge from local interaction rules. PLOS Computational Biology 18(3): e1009877. [doi.org/10.1371/journal.pcbi.1009877](https://doi.org/10.1371/journal.pcbi.1009877)
