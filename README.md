# SNO Analysis

## Conda environment
Install `conda` from script found here:
https://docs.conda.io/en/latest/miniconda.html

To create the kernel for this project:
```bash
conda env create -f environment.yml
```
To update the environment:
```bash
conda env update -f environment.yml
```

We are using `conda-lock` to ensure a reproducible environment. Please update the conda-lock files when updating the environment.

To create an environment from the lock file
```bash
conda create -n sno --file environment/conda-linux-64.lock
```

To add user to the project directory on Cheyenne:
```bash
setfacl --modify default:user:USERNAME:rwx ${project_dir}
```
