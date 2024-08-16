import os
import glob
import argparse
from ase.io import read
from host_guest.io import filetyper
from host_guest.energy import docker


def extract_energy_molecules_from_folder(host_folder, monomer_folder, number_of_host, number_of_monomers, number_of_complexes, results_folder):
    """
    A function to compute the host guest energy for host and monomer from
    their folders
    parameter
    ----------
    host_folder: str
    monomer_folder: str
    """
    seen = []
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
        json_mol_filename = os.path.join(results_folder, 'complexes.json')
        json_energy_filename = os.path.join(results_folder, 'energy.json')
    else:
        json_mol_filename = filetyper.load_data(
            os.path.join(results_folder, 'complexes.json'))
        json_energy_filename = filetyper.load_data(
            os.path.join(results_folder, 'energy.json'))
        seen = list(json_energy_filename.keys())
    new_mol = {}
    new_energy = {}
    for host_system_file in host_folder:
        for monomer_file in monomer_folder:
            monomer = read(monomer_file)
            host_system = read(host_system_file)
            host_base_name = os.path.basename(host_system_file).split('.')[0].split('_')[0]
            monomer_base_name = os.path.basename(monomer_file).split('.')[0]
            base_name = host_base_name + '_' + monomer_base_name
            if base_name not in seen:
                energy_dict, complex_molecules = docker. Dock(
                    host_system, monomer, number_of_host, number_of_monomers, number_of_complexes)
                new_mol[base_name] = complex_molecules
                new_energy[base_name] = energy_dict
                filetyper.append_json(new_energy, json_energy_filename)
                filetyper.append_json_atom(new_mol, json_mol_filename)
    return



def extract_energy_molecules_from_file(host_system_file, monomer_file, number_of_host, number_of_monomers, number_of_complexes, results_folder):
    """
    A fucntion to compute binding energy directly from files.
    parameter
    ----------
    host_folder: str
    monomer_folder: str
    """
    seen = []
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
        json_mol_filename = os.path.join(results_folder, 'complexes.json')
        json_energy_filename = os.path.join(results_folder, 'energy.json')
    else:
        json_mol_filename = filetyper.load_data(os.path.join(results_folder, 'complexes.json'))
        json_energy_filename = filetyper.load_data(os.path.join(results_folder, 'energy.json'))
        seen = list(json_energy_filename.keys())

    new_mol = {}
    new_energy = {}
    monomer = read(monomer_file)
    host_system = read(host_system_file)
    host_base_name = os.path.basename(host_system_file).split('.')[0]
    monomer_base_name = os.path.basename(monomer_file).split('.')[0]
    base_name = host_base_name + '_' + monomer_base_name
    if base_name not in seen:
        energy_dict, complex_molecules = docker. Dock(host_system, monomer, number_of_host, number_of_monomers, number_of_complexes)
        new_mol[base_name] = complex_molecules
        new_energy[base_name] = energy_dict
        filetyper.append_json(new_energy, json_energy_filename)
        filetyper.append_json_atom(new_mol, json_mol_filename)
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
