from host_guest.io import coords_library, filetyper
from host_guest.energy import docker

mof  = coords_library.load_data_as_ase('test_data/EDUSIF.cif')
molecule = coords_library.load_data_as_ase('test_data/biphenyl.xyz')

def test_docker (mof, molecule):
    energy_dict, complex_molecules = docker.Dock(mof, molecule, 1, 1, 1)
    complex_atom = complex_molecules['complex_0']
    complex_atom.write('test.cif')

test_docker (mof, molecule)