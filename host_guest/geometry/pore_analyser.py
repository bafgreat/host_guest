from scipy.optimize import brute, fmin, minimize
from ase import geometry
import numpy as np
import itertools
from ase import Atoms
from scipy.spatial import Delaunay
from host_guest.geometry import util
from scipy.spatial import ConvexHull


bondRadii, mmbond = util.BondRadii()


def set_up_for_pore(new_atom):
    '''
    Function that set up the system for pore calculation.
    This starts from the center of mass and looks for the closest atoms to the center of mass
    and computes the distance

    **parameter**
        new_atom: ase_atom

    **returns**
        COM: numpy.ndarray
        coordinates: numpy.ndarray
        elements: list
        cell: list
        pbc: bool
    '''
    com = new_atom.get_center_of_mass()
    coordinates = new_atom.positions
    elements = [i.symbol for i in new_atom]
    cell = new_atom.cell.tolist()
    pbc = False
    if len(cell) > 0:
        pbc = True

    return com, coordinates, elements, cell, pbc


def pore_diameter(com, coordinates, elements, cell, pbc):
    '''
    Function that compute the pore diameter of a system.
    This starts from the center of mass and looks for the closing atoms to the center of mass
    and computes the distance.

    **parameter**
        COM: numpy.ndarray
        coordinates: numpy.ndarray
        elements: list
        cell: list
        pbc: bool

    **returns**
        pore_d: float pore diameter
    '''
    distances = geometry.get_distances(com.reshape(1, -1), coordinates, cell, pbc=pbc)[1]

    index_closest_atom = np.argmin(distances)

    vdw_radii = bondRadii[elements[index_closest_atom]][0]

    pore_d = (distances[0][index_closest_atom] - vdw_radii)*2
    return pore_d


def correct_pore_diameter(com, *params):
    """Return negative of a pore diameter. (optimisation function)."""
    coordinates, elements, cell, pbc,  = params
    return -pore_diameter(com, coordinates, elements, cell, pbc)

def opt_pore_diameter(com, coordinates, elements, cell, pbc):
    '''
    Script that compute pore diameter by searching for optimised COM
    This is computed for assymetric systems where in the pore is not
    directly at the center of mas.
    '''

    pore_r = pore_diameter(com, coordinates, elements, cell, pbc)/2.0

    bounds = (
            (com[0]-pore_r, com[0]+pore_r),
            (com[1]-pore_r, com[1]+pore_r),
            (com[2]-pore_r, com[2]+pore_r)
        )
    minimisation = minimize(
        correct_pore_diameter, x0=com, args=(coordinates, elements, cell, pbc), bounds=bounds)

    new_com = minimisation.x

    pore_d = pore_diameter(new_com, coordinates, elements, cell, pbc)
    return pore_d, new_com


def optimise_z(z, *args):
    """Return pore diameter for coordinates optimisation in z direction."""
    x, y, coordinates, elements, cell, pbc = args

    window_com = np.array([x, y, z[0]])
    return pore_diameter(window_com, coordinates, elements, cell, pbc)


def optimise_xy(xy, *args):
    """Return negative pore diameter for x and y coordinates optimisation."""
    z, coordinates, elements, cell, pbc = args
    window_com = np.array([xy[0], xy[1], z])
    return -pore_diameter(window_com, coordinates, elements, cell, pbc)


def molecule_diameter(new_atom):
    '''
    Function that compute the maximum dimension of a molecule.
    This starts from the center of mass and looks for the closest atoms
    to the center of mass and computes the distance.

    **parameter**
        new_atom: ase_atom

    **returns**
        maxdim: float maximum dimension of the system
    '''
    dist_matrix = new_atom.get_all_distances(mic=True)
    final_matrix = np.triu(dist_matrix)
    i1, i2 = np.unravel_index(final_matrix.argmax(), final_matrix.shape)
    vdw_1 = bondRadii[new_atom[i1].symbol][0]
    vdw_2 = bondRadii[new_atom[i2].symbol][0]
    maxdim = final_matrix[i1, i2]+vdw_1 + vdw_2
    return maxdim


def pore_diameter_of_structure(ase_atom):
    '''
    Function that compute the pore diameter of a structure.

    **parameter**
        ase_atom: ase_atom
    **returns**
        pore_d: float pore diameter
    '''
    com, coordinates, elements, cell, pbc = set_up_for_pore(ase_atom)

    return opt_pore_diameter(com, coordinates, elements, cell, pbc)[0]


def calculate_euler_characteristic(atoms: Atoms):
    positions = atoms.get_positions()
    # Compute the Delaunay triangulation of the points
    delaunay = Delaunay(positions)

    # Vertices (unique points in Delaunay triangulation)
    v_ = len(delaunay.points)

    # Edges (Delaunay simplices)
    e_ = len(delaunay.convex_hull)

    # Faces (number of faces from the Delaunay simplices)
    f_ = len(delaunay.simplices)

    # Compute the Euler characteristic
    chi = v_ - e_ + f_

    return chi


def alpha_shape(points, alpha):
    """
    Compute the alpha shape (concave hull) of a set of points.
    **parameters:***
        points (ndarray): array of points.
        alpha (float): alpha value to influence the shape.
    **returns:**
        edges, triangles, volume: edges, triangles, and volume of the shape.
    """
    tetrahedra = Delaunay(points).simplices
    a, b, c, d =\
        points[tetrahedra[:, 0]], points[tetrahedra[:, 1]], points[tetrahedra[:, 2]], points[tetrahedra[:, 3]]
    normals = np.cross(b - a, c - a)
    tetrahedron_volumes = np.einsum('ij,ij->i', normals, d - a) / 6

    radius_circum_sphere =\
        np.linalg.norm(np.cross(b - a, c - a), axis=1) / (2 * np.abs(tetrahedron_volumes))
    alpha_shape_mask = radius_circum_sphere < 1 / alpha
    alpha_shape_tetrahedra = tetrahedra[alpha_shape_mask]

    edges = set()
    triangles = set()
    for simplex in alpha_shape_tetrahedra:
        for i, j in itertools.combinations(simplex, 2):
            edge = tuple(sorted((i, j)))
            edges.add(edge)
        for i, j, k in itertools.combinations(simplex, 3):
            triangle = tuple(sorted((i, j, k)))
            triangles.add(triangle)

    return edges, triangles, len(alpha_shape_tetrahedra)


def calculate_euler_characteristic(edges, triangles, num_tetrahedra):
    V = len(set(itertools.chain(*edges)))  # unique vertices
    E = len(edges)  # unique edges
    F = len(triangles)  # unique faces

    chi = V - E + F - num_tetrahedra  # Euler characteristic with volume (tetrahedra) considered
    return chi

def calculate_number_of_pores(atoms: Atoms, alpha=1.0):
    positions = atoms.get_positions()

    # Compute the alpha shape of the atomic positions
    edges, triangles, num_tetrahedra = alpha_shape(positions, alpha)

    # Compute the Euler characteristic
    chi = calculate_euler_characteristic(edges, triangles, num_tetrahedra)

    # Compute the number of pores using the relationship P = 1 - χ
    P = 1 - chi
    return P


def calculate_number_of_pores2(atoms: Atoms):
    chi = calculate_euler_characteristic(atoms)
    # Compute the number of pores using the relationship P = 1 - χ
    P = 1 - chi
    return P


def radius_from_convexhall(coords:np.ndarray):
    """
    Calculate the radius of the convex hull from a set of coordinates.


    **parameters**:
        coords (numpy.ndarray): array of coordinates.

    **returns**:
        float: radius of the convex hull.
    """
    hull = ConvexHull(coords)
    center = np.mean(coords, axis=0)
    vertices = coords[hull.vertices]
    distances = np.linalg.norm(vertices - center, axis=1)
    radius = np.min(distances)
    return radius

# Example usage with an ASE Atoms object
# atoms = Atoms(...)  # Define your Atoms object here
# num_pores = calculate_number_of_pores(atoms, alpha=1.0)
# print("Number of pores:", num_pores)


# Example usage with an ASE Atoms object
# atoms = Atoms(...)  # Define your Atoms object here
# num_pores = calculate_number_of_pores(atoms)
# print("Number of pores:", num_pores)
