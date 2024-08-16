
# host_guest
host_guest is a python module for computing binding energies between a host system and a guest molecule. It takes a host system and interpolate a series of guest molecules at random positions within a specified raddi within the host system.

For each host guest configuration, a GFN-xTB single point energy is computed and the binding energy is then derived from this computation. This is essentially a workflow to expedite the computation of binding energies.

The each complex and binding energies are then dumped in json files.


# Installation
The intsallation of host_guest requires the mamba or conda environment.

# Conda Installation
```
git clone https://github.com/bafgreat/host_guest.git
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



