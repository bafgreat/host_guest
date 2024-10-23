from host_guest.io import coords_library
from host_guest.geometry import pore_analyser

mof  = coords_library.load_data_as_ase('test_data/EDUSIF.cif')
molecule = coords_library.load_data_as_ase('test_data/biphenyl.xyz')

def test_pore_diameter(mof, molecule):
    """
    Test function to calculate the pore diameter
    """
    mole_max_diameter = pore_analyser.molecule_diameter(molecule)
    assert mole_max_diameter == 9.981643702897477
    mof_max_diameter = pore_analyser.pore_diameter_of_structure(mof)
    assert mof_max_diameter == 13.821456707655694

    print ( pore_analyser.calculate_number_of_pores(mof))

test_pore_diameter(mof, molecule)