# NICS - Nucleus Independent Chemical Shift utility: Python3, rdkit, openbabel

## Description

NICS is a program that searches for rings in input structures and places probes 1 angstrom above the center of the rings found in the structure. It also accepts manual ring input using the atom positions.

## Installation

### Step 1: Install Miniconda

Below are the command line instructions for macOS, but instructions for other operating systems are available on the [Anaconda website](https://docs.anaconda.com/miniconda/miniconda-install/).

```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
source ~/.bash_profile
conda init
```

### Step 2: Create a Virtual Environment and Activate It

It is recommended to use a virtual environment to avoid conflicting packages.

```bash
conda create -n chem-env python=3.8
conda activate chem-env
```

### Step 3: Install OpenBabel, RDKit, and Pybel Packages

Install the necessary dependencies for this program.

```bash
conda install -c conda-forge openbabel
conda install -c conda-forge rdkit
pip install pybel scipy
```

## Usage

### How to Run This Script

Make sure there is a Gaussian input file and the python script are in the directory and type:

```bash
python nics.py -I test.com
```

### Argument Description

- `-I` or `--input`: Input file name. It can be either a Gaussian16 input (_.com) or output (_.log).
- `-R` or `--rings`: Ring element indices starting from 1. Make sure single quotes are included. Example: `'1 3 6 7 9'`
- `-M` or `--method`: Method. Example: `mpwpw91`
- `-B` or `--basis`: Basis set. Example: `6-311+G(2d,p)`
- `-C` or `--charge`: Charge. Example: `-1`
- `-S` or `--spin`: Spin. Example: `2`
