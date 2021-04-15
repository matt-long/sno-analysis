# SNO Analysis

## Conda environment
Install `conda` from script found here:
https://docs.conda.io/en/latest/miniconda.html

We are using `conda-lock` to ensure a reproducible environment. Please update the conda-lock files when updating the environment.

To create an environment from the lock file
```bash
conda create -n sno --file environment/conda-linux-64.lock
```
