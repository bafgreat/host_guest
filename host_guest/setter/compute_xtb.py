import os
import subprocess
import shutil
import random
import argparse
from host_guest.io import coords_library, filetyper
from host_guest.energy.compute_sp import compute_xtb_energy, compute_energy_of_atom

def run_xtb_energy(input_system, output_folder):
    """
    A function to run xtb energy calculation on a system.

    Parameters
    ----------
    input_system: file or directory
        Input system to be computed.
    output_folder: str
        Output folder for storing the results.
    """

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if os.path.isdir(input_system):
        compute_energy_of_atom(input_system, output_folder)
    else:
        ase_atoms = coords_library.read_and_return_ase_atoms(input_system)
        base_name = os.path.basename(input_system).split('.')[0]
        energy = compute_xtb_energy(ase_atoms)
        filetyper.write_json(energy, f'{output_folder}/{base_name}_energy.json')


def main():
    '''
    Command line interface for creating cifs from
    '''
    parser = argparse.ArgumentParser(
        description="""
        CLI for computing the xtb energy of a system from folder or file.
        if the system is a file, it reads and calculates the energy directly
        if the system is a directory, it loops through all files and calculates
        the energy for each file.
        """
        )
    parser.add_argument('input_system', type=str,
                        help="""system containg coordinates. It could be folder containing
                        structure file or it could be a file""")

    parser.add_argument('-r', '--results_folder', type=str,
                        default='energy_folders', help='directory to save cif files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print verbose output')
    args = parser.parse_args()
    run_xtb_energy(args.input_system, args.results_folder)
