from host_guest.io import coords_library, filetyper
from host_guest.energy import docker

mof  = coords_library.load_data_as_ase('test_data/First_MOF.xyz')
molecule = coords_library.load_data_as_ase('test_data/biphenyl.xyz')

def test_docker (mof, molecule):
    energy_dict, complex_molecules = docker.Dock(mof, molecule, 1, 1, 5)
    for i, complexes in enumerate(complex_molecules):
        complex_atom = complex_molecules[complexes]
        complex_atom.write(f'test_{i}.mol')

test_docker (mof, molecule)