from scipy.optimize import brute, fmin, minimize
from ase import geometry
import numpy as np
from host_guest.geometry import util

bondRadii, mmbond = util.BondRadii()

def set_up_for_pore(new_atom):
    '''
    Function that set up the system for pore calculation.
    This starts from the center of mass and looks for the closest atoms to the center of mass
    and computes the distance
    parameter
    ---------
    new_atom: ase_atom
    Returns
    COM: numpy.ndarray
    coordinates: numpy.ndarray
    elements: list
    cell: list
    pbc: bool
    '''
    com = new_atom.get_center_of_mass()
    coordinates =new_atom.positions
    elements = [i.symbol for i in new_atom]
    cell =  new_atom.cell.tolist()
    pbc = False
    if len(cell)>0:
        pbc=True

    return com, coordinates, elements, cell, pbc

def pore_diameter(COM, coordinates, elements, cell, pbc):
    '''
    Function that compute the pore diameter of a system.
    This starts from the center of mass and looks for the closing atoms to the center of mass
    and computes the distance.
    parameter:
    COM: numpy.ndarray
    coordinates: numpy.ndarray
    elements: list
    cell: list
    pbc: bool
    Returns
    pore_d: float pore diameter
    '''
    distances = geometry.get_distances(COM.reshape(1, -1), coordinates, cell, pbc=pbc)[1]


    index_closest_atom = np.argmin(distances)

    vdw_radii = bondRadii[ elements[index_closest_atom]][0]

    pore_d = (distances[0][index_closest_atom ]- vdw_radii )*2
    return  pore_d

def correct_pore_diameter(COM, *params):
    """Return negative of a pore diameter. (optimisation function)."""
    coordinates, elements, cell, pbc,  = params
    return -pore_diameter(COM, coordinates, elements, cell, pbc)

def opt_pore_diameter(COM, coordinates, elements, cell, pbc):
    '''
    Script that compute pore diameter by searching for optimised COM
    This is computed for assymetric systems where in the pore is not
    directly at the center of mas.
    '''

    pore_r = pore_diameter(COM, coordinates, elements, cell, pbc)/2.0

    bounds = (
            (COM[0]-pore_r, COM [0]+pore_r),
            (COM[1]-pore_r, COM[1]+pore_r),
            (COM[2]-pore_r, COM [2]+pore_r)
        )
    minimisation =  minimize(
        correct_pore_diameter, x0=COM, args=(coordinates,elements, cell, pbc ), bounds=bounds)

    new_COM = minimisation.x

    pore_d = pore_diameter(new_COM, coordinates, elements, cell, pbc)
    return pore_d,  new_COM

def optimise_z(z, *args):
    """Return pore diameter for coordinates optimisation in z direction."""
    x, y, coordinates, elements, cell, pbc= args

    window_com = np.array([x, y, z[0]])
    return pore_diameter(window_com, coordinates, elements, cell, pbc)

def optimise_xy(xy, *args):
    """Return negative pore diameter for x and y coordinates optimisation."""
    z, coordinates, elements, cell, pbc= args
    window_com = np.array([xy[0], xy[1], z])
    return -pore_diameter(window_com, coordinates, elements, cell, pbc)

def molecule_diameter(new_atom):
    '''
    Function that compute the maximum dimension of a molecule.
    This starts from the center of mass and looks for the closest atoms to the center of mass
    and computes the distance.
    parameter:
    new_atom: ase_atom
    Returns
    maxdim: float maximum dimension of the system
    '''
    dist_matrix = new_atom.get_all_distances(mic=True)
    final_matrix = np.triu(dist_matrix)
    i1, i2 = np.unravel_index(final_matrix.argmax(), final_matrix.shape)
    vdw_1 = bondRadii[ new_atom[i1].symbol][0]
    vdw_2 = bondRadii[ new_atom[i2].symbol][0]
    maxdim = final_matrix[i1, i2]+vdw_1 +vdw_2
    return maxdim

def pore_diameter_of_structure(ase_atom):
    '''
    Function that compute the pore diameter of a structure.
    parameter:
    ase_atom: ase_atom
    Returns
    pore_d: float pore diameter
    '''
    com, coordinates, elements, cell, pbc = set_up_for_pore(ase_atom)

    return opt_pore_diameter(com, coordinates, elements, cell, pbc)[0]