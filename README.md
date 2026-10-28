# CannTeen

Analysis pipeline for the **CannTeen** study, investigating MRS (magnetic resonance spectroscopy) outcomes in adolescents and adults under acute THC and THC+CBD administration.

## Pipeline

This repository contains the main statistical analysis script:

- **[CannTeen_Analysis.R](CannTeen_Analysis.R)** — R script that merges clinical, questionnaire, and drug-condition data and runs the main statistical analyses (linear mixed-effects models, pairwise/permutation tests, correlations).

The upstream imaging preprocessing steps (MRS voxel segmentation in MATLAB and a voxel-position/density check in Python) are maintained locally and are not published here.

## Environment

`CannTeen_Analysis.R` requires R with the packages listed at the top of the file (e.g. `lmerTest`, `emmeans`, `ptestr`).

[Voxel_Check_env.yml](Voxel_Check_env.yml) is a conda environment file kept for reference; it supports the local (unpublished) Python voxel-position/density check step and is not required to run `CannTeen_Analysis.R`.

## Data

Raw and participant-level data are not included in this repository and must be supplied separately; file paths in the script will need to be updated to match your local setup.

## License

This project is licensed under the [MIT License](LICENSE).

