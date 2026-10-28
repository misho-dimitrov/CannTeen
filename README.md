# CannTeen

Analysis pipeline for the **CannTeen** study, investigating MRS (magnetic resonance spectroscopy) outcomes in adolescents and adults under acute THC and THC+CBD administration.

## Pipeline

Scripts are numbered in the order they are intended to be run:

1. **[1_MRS_Segmentation.m](1_MRS_Segmentation.m)** — MATLAB script that uses [Gannet](https://github.com/markmikkelsen/Gannet) and FSL to segment MRS voxels into grey matter, white matter, and CSF for each subject/session.
2. **[2_MRS_Voxel_Position_Check.py](2_MRS_Voxel_Position_Check.py)** — Python script that loads per-subject voxel masks, computes group-level voxel density overlap maps, and plots them on a glass brain.
3. **[3_Analysis.R](3_Analysis.R)** — R script that merges clinical, questionnaire, and drug-condition data and runs the main statistical analyses (linear mixed-effects models, pairwise/permutation tests, correlations).

## Environment

The Python dependencies for script 2 are listed in [Voxel_Check_env.yml](Voxel_Check_env.yml). Create the conda environment with:

```bash
conda env create -f Voxel_Check_env.yml
conda activate Voxel_Check_env
```

Script 1 additionally requires MATLAB with [FSL](https://fsl.fmrib.ox.ac.uk/) and [Gannet](https://github.com/markmikkelsen/Gannet) on the path. Script 3 requires R with the packages listed at the top of the file (e.g. `lmerTest`, `emmeans`, `ptestr`).

## Data

Raw and participant-level data are not included in this repository and must be supplied separately; file paths in the scripts will need to be updated to match your local setup.

## License

This project is licensed under the [MIT License](LICENSE).
