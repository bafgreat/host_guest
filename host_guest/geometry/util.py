
from ase.data import chemical_symbols,covalent_radii
import numpy as np
from copy import deepcopy
from sklearn.metrics.pairwise import euclidean_distances
from ase import geometry

def BondRadii():
    bondRadii = {symbol:[r, r-0.1, r-0.2] for symbol,r in zip(chemical_symbols,covalent_radii)}
    mmbond = {'Cr':4,'Mo':4,'V':3,'Rh':1,'Ru':2,'W':4,'Os':3,'Re':4,'Pt':1,'Tc':4,'Ir':1,'Pd':1}
    return bondRadii, mmbond


def sphere_volume(sphere_radius):
    """Return volume of a sphere."""
    return (4 / 3 * np.pi * sphere_radius**3)


def calc_triangle_normal(p1, p2, p3):
    v1 = p2 - p1
    v2 = p3 - p1
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    d = np.dot(v1, v2)
    c = np.cross(v1, v2)

    return c

def asphericity(S):
    return (S[0] - (S[1] + S[2]) / 2)


def acylidricity(S):
    return (S[1] - S[2])


def relative_shape_anisotropy(S):
    return (1 - 3 * (
        (S[0] * S[1] + S[0] * S[2] + S[1] * S[2]) / (np.sum(S))**2))



def get_tensor_eigenvalues(T, sort=False):
    if sort:
        return (sorted(np.linalg.eigvals(T), reverse=True))
    else:
        return (np.linalg.eigvals(T))

def get_gyration_tensor(new_atom):
    """
    Return the gyration tensor of a molecule.
    The gyration tensor should be invariant to the molecule's position.
    The known formulas for the gyration tensor have the correction for the
    centre of mass of the molecule, therefore, the coordinates are first
    corrected for the centre of mass and essentially shifted to the origin.
    Parameters
    ----------
    new_atom : ASE atom type
    Returns
    -------
    numpy.ndarray
        The gyration tensor of a molecule invariant to the molecule's position.
    """
    #Compute centre of mass
    COM = new_atom.get_center_of_mass()
    #Move all coordinates to the centre of mass
    coordinates = new_atom.positions - COM
   # Calculate diagonal and then other values of the matrix.
    diag = np.sum(coordinates**2, axis=0)
    xy = np.sum(coordinates[:, 0] * coordinates[:, 1])
    xz = np.sum(coordinates[:, 0] * coordinates[:, 2])
    yz = np.sum(coordinates[:, 1] * coordinates[:, 2])
    S = np.array([[diag[0], xy, xz], [xy, diag[1], yz],
                  [xz, yz, diag[2]]]) / coordinates.shape[0]
    return (S)


def get_inertia_tensor(new_atom):
    """
    Return the tensor of inertia a molecule.
    Parameters
    ----------
    new_atom: ASE atom types
    Returns
    -------
    numpy.ndarray
        The tensor of inertia of a molecule.
    """
    coordinates = new_atom.positions
    pow2 = coordinates**2
    molecular_weight = new_atom.get_masses()


    diag_1 = np.sum(molecular_weight * (pow2[:, 1] + pow2[:, 2]))
    diag_2 = np.sum(molecular_weight * (pow2[:, 0] + pow2[:, 2]))
    diag_3 = np.sum(molecular_weight * (pow2[:, 0] + pow2[:, 1]))

    mxy = np.sum(-molecular_weight * coordinates[:, 0] * coordinates[:, 1])
    mxz = np.sum(-molecular_weight * coordinates[:, 0] * coordinates[:, 2])
    myz = np.sum(-molecular_weight * coordinates[:, 1] * coordinates[:, 2])

    inertia_tensor = np.array([[diag_1, mxy, mxz], [mxy, diag_2, myz],
                               [mxz, myz, diag_3]]) / coordinates.shape[0]
    print (inertia_tensor)

    return inertia_tensor

def normalize_vector(vector):
    """
    Normalize a vector.
    A new vector is returned, the original vector is not modified.
    Parameters
    ----------
    vector : np.array
        The vector to be normalized.
    Returns
    -------
    np.array
        The normalized vector.
    """
    v = np.divide(vector, np.linalg.norm(vector))
    return np.round(v, decimals=4)

def rotation_matrix_arbitrary_axis(angle, axis):
    """
    Return a rotation matrix of `angle` radians about `axis`.
    Parameters
    ----------
    angle : int or float
        The size of the rotation in radians.
    axis : numpy.array
        A 3 element aray which represents a vector. The vector is the
        axis about which the rotation is carried out.
    Returns
    -------
    numpy.array
        A 3x3 array representing a rotation matrix.
    """
    axis = normalize_vector(axis)

    a = np.cos(angle / 2)
    b, c, d = axis * np.sin(angle / 2)

    e11 = np.square(a) + np.square(b) - np.square(c) - np.square(d)
    e12 = 2 * (b * c - a * d)
    e13 = 2 * (b * d + a * c)

    e21 = 2 * (b * c + a * d)
    e22 = np.square(a) + np.square(c) - np.square(b) - np.square(d)
    e23 = 2 * (c * d - a * b)

    e31 = 2 * (b * d - a * c)
    e32 = 2 * (c * d + a * b)
    e33 = np.square(a) + np.square(d) - np.square(b) - np.square(c)

    return np.array([[e11, e12, e13], [e21, e22, e23], [e31, e32, e33]])


def principal_axes(new_atom):
    return (np.linalg.eig(get_inertia_tensor(new_atom))[1].T)

def align_principal_ax(new_atom):
    """
    Align the principal axes of a molecule to the Cartesian axes.
    Parameters
    ----------
    new_atom : ASE Atom object
        The molecule to be aligned.
    Returns
    -------
    (numpy.ndarray, list)
        The aligned coordinates and rotation matrices.
    """
    coor = deepcopy(new_atom.positions)
    new_coor = []
    rot = []
    elements = np.array([i.symbol for  i in new_atom])
    print (elements)
    for i, j in zip([2, 1, 0], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]):
        p_axes = principal_axes(new_atom)
        r_vec = np.cross(p_axes[i], np.array(j))
        sin = np.linalg.norm(r_vec)
        cos = np.dot(p_axes[i], np.array(j))
        ang = np.arctan2(sin, cos)
        R_mat = np.matrix(rotation_matrix_arbitrary_axis(ang, r_vec))
        rot.append(R_mat)
        for i in coor:
            new_coord = R_mat * i.reshape(-1, 1)
            new_coor.append(np.array(new_coord.reshape(1, -1))[0])
        new_coor = np.array(new_coor)
        coor = new_coor
        new_coor = []
    return (coor, rot)

def vector_analysis(vector, coordinates, elements_vdw, cell , pbc=True, increment=1.0):
    """Analyse a sampling vector's path for window analysis purpose."""
    # Calculate number of chunks if vector length is divided by increment.
    chunks = int(np.linalg.norm(vector) // increment)
    #print (chunks)
    # Create a single chunk.
    chunk = vector / chunks
    # Calculate set of points on vector's path every increment.
    vector_pathway = np.array([chunk * i for i in range(chunks + 1)])
    #test = vector_pathway[0].reshape(1, -1)
    #test2 = geometry.get_distances(test,coordinates, cell , pbc=True)
    analysed_vector = np.array([
        np.amin(
            geometry.get_distances(i.reshape(1, -1), coordinates, cell , pbc=True)[1] - elements_vdw)
        for i in vector_pathway
    ])
    if all(i > 0 for i in analysed_vector):
        pos = np.argmin(analysed_vector)
        # As first argument we need to give the distance from the origin.
        dist = np.linalg.norm(chunk * pos)
        return np.array([dist, analysed_vector[pos] * 2, *chunk * pos, *vector])


def vector_preanalysis(vector, coordinates, elements_vdw,cell, pbc=True, increment=1.0):
    norm_vec = vector/np.linalg.norm(vector)
    intersections = []
    origin = np.mean(coordinates, axis=0)
    L = coordinates - origin
    t_ca = np.dot(L, norm_vec)
    d = np.sqrt(np.einsum('ij,ij->i', L, L) - t_ca**2)

    under_sqrt = elements_vdw**2 - d**2

    previous = under_sqrt .reshape(1,-1)
    diag = np.diagonal((previous)*(previous).T)

    #diag = under_sqrt.diagonal()
    positions = np.argwhere(diag > 0)
    for pos in positions:
        t_hc = np.sqrt(diag[pos[0]])
        t_0 = t_ca[pos][0] - t_hc
        t_1 = t_ca[pos][0] + t_hc

        P_0 = origin + np.dot(t_0, norm_vec)
        P_1 = origin + np.dot(t_1, norm_vec)
        #print(np.linalg.norm(P_0), np.linalg.norm(P_1))
        if np.linalg.norm(P_0) < np.linalg.norm(P_1):
            intersections.append(1)
        else:
            intersections.append(0)
    if sum(intersections) == 0:
        return vector_analysis(vector, coordinates, elements_vdw, cell , pbc=True, increment=1.0)

def normal_vector(origin, vectors):
    """Return normal vector for two vectors with same origin."""
    return np.cross(vectors[0] - origin, vectors[1] - origin)

def angle_between_vectors(x, y):
    """Calculate the angle between two vectors x and y."""
    first_step = abs(x[0] * y[0] + x[1] * y[1] + x[2] * y[2]) / (
        np.sqrt(x[0]**2 + x[1]**2 + x[2]**2) *
        np.sqrt(y[0]**2 + y[1]**2 + y[2]**2))
    second_step = np.arccos(first_step)
    return (second_step)

def shift_com(COM, coordinates, elements,  com_adjust=np.zeros(3)):
    """
    Return coordinates translated by some vector.

    Parameters
    ----------
    elements : numpy.ndarray
        An array of all elements (type: str) in a molecule.

    coordinates : numpy.ndarray
        An array containing molecule's coordinates.

    com_adjust : numpy.ndarray (default = [0, 0, 0])

    Returns
    -------
    numpy.ndarray
        Translated array of molecule's coordinates.

    """

    com = np.array([COM - com_adjust] * coordinates.shape[0])
    return coordinates - com


def vector_analysis_reversed(vector, coordinates, elements_vdw):
    norm_vec = vector/np.linalg.norm(vector)
    intersections = []
    origin = np.mean(coordinates, axis=0)
    L = coordinates - origin
    t_ca = np.dot(L, norm_vec)
    d = np.sqrt(np.einsum('ij,ij->i', L, L) - t_ca**2)
    under_sqrt = (elements_vdw**2 - d**2).reshape(1,-1)

    diag = under_sqrt.diagonal()
    positions = np.argwhere(diag > 0)
    for pos in positions:
        t_hc = np.sqrt(diag[pos[0]])
        t_0 = t_ca[pos][0] - t_hc
        t_1 = t_ca[pos][0] + t_hc

        P_0 = origin + np.dot(t_0, norm_vec)
        P_1 = origin + np.dot(t_1, norm_vec)
        if np.linalg.norm(P_0) < np.linalg.norm(P_1):
            intersections.append([np.linalg.norm(P_1),  P_1])
    if intersections:
        intersection = sorted(intersections, reverse=True)[0][1]
        dist_origin = np.linalg.norm(intersection)
        return [dist_origin, intersection]


def vector_analysis_pore_shape(vector, coordinates, elements_vdw):
    norm_vec = vector/np.linalg.norm(vector)
    intersections = []
    origin = np.mean(coordinates, axis=0)
    L = coordinates - origin
    t_ca = np.dot(L, norm_vec)
    d = np.sqrt(np.einsum('ij,ij->i', L, L) - t_ca**2)
    under_sqrt = (elements_vdw**2 - d**2).reshape(1,-1)
    diag = under_sqrt.diagonal()
    positions = np.argwhere(diag > 0)
    for pos in positions:
        t_hc = np.sqrt(diag[pos[0]])
        t_0 = t_ca[pos][0] - t_hc
        t_1 = t_ca[pos][0] + t_hc

        P_0 = origin + np.dot(t_0, norm_vec)
        P_1 = origin + np.dot(t_1, norm_vec)
        # print(np.linalg.norm(P_0), np.linalg.norm(P_1))
        if np.linalg.norm(P_0) < np.linalg.norm(P_1):
            intersections.append([np.linalg.norm(P_0), P_0])
    if intersections:
        return sorted(intersections)[0][1]
