import os
import glob
from ase.io import read
from ase import Atoms
import argparse
from host_guest.io import filetyper, coords_library
from host_guest.energy import docker

def get_ase_atoms(js_data):
    """
    Convert a JSON data containing complexes to ASE atoms.

    Parameters
    ----------
    js_data: dict
        Dictionary containing complexes.

    Returns
    -------
    ase_atoms: ase.Atoms
        ASE atoms object representing the complexes.
    """
    ase_atoms = Atoms(symbols=js_data.get("labels"), positions=js_data.get('positions'), cell=js_data.get('lattice_vectors'), pbc=js_data.get('periodic'))
    return ase_atoms

def json_to_cif(json_file, output_dir):
    """
    Convert a JSON file containing complexes to CIF files.

    Parameters
    ----------
    json_file: str
        Path to the JSON file containing complexes.
    output_dir: str
        Path to the output directory where CIF files will be stored.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    json_data = filetyper.load_data(json_file)

    for dirname in json_data:
        path_to_save = os.path.join(output_dir, dirname)
        if not os.path.exists(path_to_save):
            os.makedirs(path_to_save)
        for base_name in json_data.get(dirname):
            cif_file = os.path.join(path_to_save, base_name + '.cif')
            atom_data = json_data[dirname][base_name]
            ase_atoms = get_ase_atoms(atom_data)
            ase_atoms.write(cif_file)


def main():
    '''
    Command line interface for creating cifs from
    '''
    parser = argparse.ArgumentParser(
        description='CLI for creating a cif folder from json complex data files')
    parser.add_argument('json_file', type=str,
                        help="""json file containing ase atoms or complexes.
                        NB: This should be in the same format as the want created
                        by the host_guest module. This will not work if the
                        json file is formattted differently.""")

    parser.add_argument('-r', '--results_folder', type=str,
                        default='Complex_cif_folders', help='directory to save cif files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print verbose output')
    args = parser.parse_args()

    if os.path.isfile(args.json_file):
        json_to_cif(args.json_file, args.results_folder)
    elif os.path.isdir(args.json_file):
        all_json = [os.path.join(args.json_file, f) for f in os.listdir(
            args.json_file) if f.endswith('.json')]
        for js_file in all_json:
            try:
                json_to_cif(js_file, args.results_folder)
            except Exception as e:
                print(f'Error in {js_file}: {e}')
                pass

