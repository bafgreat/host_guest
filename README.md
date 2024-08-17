![intro](images/intro.gif)

# `host_guest` Python Module

## Overview

`host_guest` is a Python module for generating complexes between a host system and a guest molecule. This module facilitates the process by automating the placement of guest molecules within a specified radius around the host system and calculating the binding energies using the GFN-xTB method.

## Features

- **Automated Host-Guest Configuration**:
  The module accepts a host system and guest molecules as input and automatically places a series of guest molecules at random positions within a radius around the host. One or more host or guest molecules can be inserted for each configuration.

- **Energy Calculation**:
  For each generated host-guest configuration, the module computes the single-point energy using the GFN-xTB method.

- **Data Output**:
  All computed binding energies and the corresponding host-guest configurations are saved in JSON files consequently providing an organized and accessible format for further analysis.

# Installation

The intsallation of host_guest requires mamba or conda environment.This is because `xtb` requires `conda-forges` which makes it work on all platforms.

# Conda Installation

This document provides a step-by-step guide for installing the `host_guest` project using Conda. Follow the instructions below to set up your environment.

## Prerequisites

Before you begin, ensure that you have the following installed on your system:

- **Conda**: You can install Anaconda or Miniconda, which provides the Conda package manager.
- **Git**: Make sure Git is installed to clone the repository.

## Installation Steps

1. **Clone the Repository**

Use the following command to clone the `host_guest` repository from GitHub:

```
git clone https://github.com/bafgreat/host_guest.git
cd host_guest
conda env create -f  conda-env.yml
```

The above line of code will create a new conda environment called
host_guest and will also install every dependencies.

After installation activate your environment as follows:

```
conda activate host_guest
```

# Usage

## Loading system

So far all ase readable files, qchem, AMS, Gaussian input and output files can be loaded directly.

```
from host_guest.io import coords_library
mof  = coords_library.load_data_as_ase('test_data/EDUSIF.cif')
molecule = coords_library.load_data_as_ase('test_data/biphenyl.xyz')
```

## Computing maximum molecule diameter

You can compute the maximum diameter of your molecule to verify whether the whether it can fit inside a pore

```
from host_guest.io import coords_library
from host_guest.geometry import pore_analyser
molecule = coords_library.load_data_as_ase('test_data/biphenyl.xyz')
mole_max_diameter = pore_analyser.molecule_diameter(molecule)
```

## Maximum diameter of pore

Compute the diameter of a porous system

```
from host_guest.io import coords_library
from host_guest.geometry import pore_analyser
mof  = coords_library.load_data_as_ase('test_data/EDUSIF.cif')
mof_max_diameter = pore_analyser.pore_diameter_of_structure(mof)
```

## Binding energy

```
from host_guest.io import coords_library
from host_guest.energy import docker

mof  = coords_library.load_data_as_ase('test_data/EDUSIF.cif')
molecule = coords_library.load_data_as_ase('test_data/biphenyl.xyz')
energy_dict, complex_molecules = docker.Dock(mof, molecule, 1, 1, 1)
print (energy_dict)
```
