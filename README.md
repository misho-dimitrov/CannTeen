# CannTeen: The acute effects of cannabis with and without cannabidiol in adults and adolescents

Analysis pipeline for an investigation of MRS (magnetic resonance spectroscopy) outcomes in adults and adolescents under acute THC and THC+CBD administration, part of the [**CannTeen** study](https://onlinelibrary.wiley.com/doi/full/10.1111/add.16154).

## Pipeline

This repository contains the main statistical analysis script:

- **[CannTeen_Analysis.R](CannTeen_Analysis.R)** — R script that merges clinical, questionnaire, and drug-condition data, cleans and QC-filters MRS sessions, then fits linear mixed-effects models of Glx by condition and age group with permutation-based significance and partial R² (plus post-hoc pairwise comparisons and estimated marginal means). It also runs correlation analyses between Glx shifts and THC/CBD absorption, and generates the accompanying plots.

## Data

Raw and participant-level data are not included in this repository and must be supplied separately; file paths in the script will need to be updated to match your local setup.

## License

This project is licensed under the [MIT License](LICENSE).
