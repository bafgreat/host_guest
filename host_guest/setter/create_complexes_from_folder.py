import os
import glob
import argparse
from ase.io import read
from host_guest.io import filetyper, coords_library
from host_guest.energy import docker


def extract_energy_molecules_from_folder(host_folder, monomer_folder, number_of_host, number_of_monomers, number_of_complexes, results_folder):
    """
    A function to compute the host-guest energy for host and monomer from
    their folders.
    Parameters
    ----------
    host_folder: list of str
        List of paths to host files.
    monomer_folder: list of str
        List of paths to monomer files.
    number_of_host: int
        Number of host molecules.
    number_of_monomers: int
        Number of monomer molecules.
    number_of_complexes: int
        Number of complexes to consider.
    results_folder: str
        Path to the folder where results are stored.
    """
    seen = []

    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    json_mol_filename = os.path.join(results_folder, 'complexes.json')
    json_energy_filename = os.path.join(results_folder, 'energy.json')

    if os.path.exists(json_mol_filename):
        json_mol_data = filetyper.load_data(json_mol_filename)
    else:
        json_mol_data = {}

    if os.path.exists(json_energy_filename):
        json_energy_data = filetyper.load_data(json_energy_filename)
        seen = list(json_energy_data.keys())
    else:
        json_energy_data = {}


    for host_system_file in host_folder:
        for monomer_file in monomer_folder:
            monomer = coords_library.load_data_as_ase(monomer_file)
            host_system = coords_library.load_data_as_ase(host_system_file)
            host_base_name = os.path.basename(host_system_file).split('.')[0].split('_')[0]
            monomer_base_name = os.path.basename(monomer_file).split('.')[0]
            base_name = host_base_name + '_' + monomer_base_name

            if base_name not in seen:
                energy_dict, complex_molecules = docker.Dock(
                    host_system, monomer, number_of_host, number_of_monomers, number_of_complexes
                )
                json_mol_data[base_name] = complex_molecules
                json_energy_data[base_name] = energy_dict

                filetyper.append_json(json_energy_data, json_energy_filename)
                filetyper.append_json_atom(json_mol_data, json_mol_filename)

    return

def extract_energy_molecules_from_file(host_system_file, monomer_file, number_of_host, number_of_monomers, number_of_complexes, results_folder):
    """
    A function to compute binding energy directly from files.
    Parameters
    ----------
    host_system_file: str
        Path to the host file.
    monomer_file: str
        Path to the monomer file.
    number_of_host: int
        Number of host molecules.
    number_of_monomers: int
        Number of monomer molecules.
    number_of_complexes: int
        Number of complexes to consider.
    results_folder: str
        Path to the folder where results are stored.
    """
    seen = []

    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    json_mol_filename = os.path.join(results_folder, 'complexes.json')
    json_energy_filename = os.path.join(results_folder, 'energy.json')

    if os.path.exists(json_mol_filename):
        json_mol_data = filetyper.load_data(json_mol_filename)
    else:
        json_mol_data = {}

    if os.path.exists(json_energy_filename):
        json_energy_data = filetyper.load_data(json_energy_filename)
        seen = list(json_energy_data.keys())
    else:
        json_energy_data = {}


    monomer = coords_library.load_data_as_ase(monomer_file)
    host_system = coords_library.load_data_as_ase(host_system_file)
    host_base_name = os.path.basename(host_system_file).split('.')[0]
    monomer_base_name = os.path.basename(monomer_file).split('.')[0]
    base_name = host_base_name + '_' + monomer_base_name

    if base_name not in seen:
        energy_dict, complex_molecules = docker.Dock(
            host_system, monomer, number_of_host, number_of_monomers, number_of_complexes
        )
        json_mol_data[base_name] = complex_molecules
        json_energy_data[base_name] = energy_dict

        filetyper.append_json(json_energy_data, json_energy_filename)
        filetyper.append_json_atom(json_mol_data, json_mol_filename)

    return



def main():
    '''
    Command line interface for computing docker
    '''
    parser = argparse.ArgumentParser(
        description='Run work_flow function with optional verbose output')
    parser.add_argument('host_folder', type=str,
                        help='A folder containing host systems')

    parser.add_argument('monomer_folder', type=str,
                        help='A folder containing guest systems')

    parser.add_argument('-nh', '--number_of_host', type=int,
                        default=1, help='The number of host systems')

    parser.add_argument('-nm', '--number_of_monomers', type=int,
                        default=1, help='The number of monomer systems')

    parser.add_argument('-nc', '--number_of_complexes', type=int,
                        default=1, help='The number of monomer systems')

    parser.add_argument('-r', '--results_folder', type=str,
                        default='results_folders', help='directory to save output files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print verbose output')
    args = parser.parse_args()

    host_folder = [os.path.join(args.host_folder, f) for f in os.listdir(
        args.host_folder) if f.endswith('.cif')]

    monomer_folder = [os.path.join(args.monomer_folder, f) for f in os.listdir(
        args.monomer_folder) if f.endswith('.xyz')]

    extract_energy_molecules_from_folder(host_folder, monomer_folder, args.number_of_host,
                                         args.number_of_monomers, args.number_of_complexes, args.results_folder)

def main2():
    '''
    Command line interface for computing docker
    '''
    parser = argparse.ArgumentParser(
        description='Run work_flow function with optional verbose output')
    parser.add_argument('host_file', type=str,
                        help='A file containing host systems')

    parser.add_argument('monomer_file', type=str,
                        help='A file containing guest systems')

    parser.add_argument('-nh', '--number_of_host', type=int,
                        default=1, help='The number of host systems')

    parser.add_argument('-nm', '--number_of_monomers', type=int,
                        default=1, help='The number of monomer systems')

    parser.add_argument('-nc', '--number_of_complexes', type=int,
                        default=1, help='The number of monomer systems')

    parser.add_argument('-r', '--results_folder', type=str,
                        default='results_folders', help='directory to save output files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print verbose output')
    args = parser.parse_args()


    extract_energy_molecules_from_file(args.host_file, args.monomer_file, args.number_of_host, args.number_of_monomers, args.number_of_complexes, args.results_folder)
