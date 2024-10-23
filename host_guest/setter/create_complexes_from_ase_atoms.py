import os
import glob
from ase.io import read
import argparse
from host_guest.io import filetyper, coords_library
from host_guest.energy import docker


def extract_energy_molecules_from_file(list_of_hosts, list_of_monomers, number_of_host, number_of_monomers, number_of_complexes, selector, results_folder):
    """
    A fucntion to extract the energy of a given host-guest system.
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
    if os.path.exists(json_mol_filename):
        json_mol_filename = filetyper.load_data(
            os.path.join(results_folder, 'complexes.json'))
        seen = list(json_energy_filename.keys())

    if os.path.exists(json_energy_filename):
        json_energy_filename = filetyper.load_data(
            os.path.join(results_folder, 'energy.json'))


    new_mol = {}
    new_energy = {}

    if selector is not None:
        list_of_hosts = sorted(glob.glob(f'{list_of_hosts}/{selector}*'))
    else:
        list_of_hosts = sorted([os.path.join(list_of_hosts, i)
                               for i in os.listdir(list_of_hosts)])
    list_of_monomers = sorted([os.path.join(list_of_monomers, i)
                              for i in os.listdir(list_of_monomers)])

    for host_system_file in list_of_hosts:
        for monomer_file in list_of_monomers:
            monomer = coords_library.load_data_as_ase(monomer_file)
            host_system = coords_library.load_data_as_ase(host_system_file)
            host_base_name = os.path.basename(host_system_file).split('.')[0]
            monomer_base_name = os.path.basename(monomer_file).split('.')[0]
            base_name = host_base_name + '_' + monomer_base_name
            if base_name not in seen:
                energy_dict, complex_molecules = docker. Dock(
                    host_system, monomer, number_of_host, number_of_monomers, number_of_complexes)
                new_mol[base_name] = complex_molecules
                new_energy[base_name] = energy_dict
                filetyper.append_json(new_energy, json_energy_filename)
                filetyper.append_json_atom(new_mol, json_mol_filename)
            else:
                print(f'{base_name} has already been computed')
    return


def main():
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

    parser.add_argument('-sl', '--selector', type=str,  default='None',
                        help='selector to select hosts')

    parser.add_argument('-r', '--results_folder', type=str,
                        default='results_folders', help='directory to save output files')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print verbose output')
    args = parser.parse_args()

    extract_energy_molecules_from_file(args.host_file, args.monomer_file, args.number_of_host,
                                       args.number_of_monomers, args.number_of_complexes, args.selector,  args.results_folder)
