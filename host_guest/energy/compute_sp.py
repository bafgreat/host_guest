import os
import subprocess
import shutil
import random
from host_guest.io import coords_library, filetyper

def compute_xtb_energy(ase_atoms):
    """
    A function that computes the xtb energy of any ase atom.
    parameter
    ----------
    ase_atoms: ase.Atoms object
    """
    base_dir = os.getcwd()
    random_number = random.uniform(0, 10**8)
    result_folder = 'tmp_dir_'+str(random_number)
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    os.chdir(result_folder)
    tmp_in = 'tmp_energy.gen'
    tmp_out = 'tmp_energy.out'
    coords_library.write_ase_atoms(ase_atoms, tmp_in)
    os.system(f'xtb --sp --gfn 2 --tblite --spinpol {tmp_in} > {tmp_out}')
    energy = read_xtb_energy(tmp_out)
    os.chdir(base_dir)
    if os.path.exists(result_folder):
        shutil.rmtree(result_folder)
    return energy


def read_xtb_energy(filename):
    """
    A function that opens the output file
    and reads the xtb energy

    **parameter**
        filename: output file

    **return**
        A dictionary containing the xtb energy
    """
    contents = filetyper.get_contents(filename)
    energy = {}
    for line in contents:
        if 'TOTAL ENERGY' in line:
            data = line.split()
            energy['energy_kcal_mol'] = float(data[3])*627.5095
        if 'HOMO-LUMO GAP' in line:
            data = line.split()
            energy['homo_lumo_ev'] = float(data[3])
    return energy


def compute_energy_of_atom(folder_of_atoms, output_folder):
    """
    A function to compute the energy of monomers to be inserted into the host.
    The function goes through a folder and searches for all monomers
    and then computes their single point energies.

    **parameter**
        folder_of_atoms: A folder containing all the monomers
        output_folder: The output folder to store the results
    """
    all_energy = {}
    json_filename = f'{output_folder}/energy_of_atoms.json'
    for filename in folder_of_atoms:
        basename = filename.split('/')[-1].split('_')[0]
        ase_atoms = coords_library.read_and_return_ase_atoms(filename)
        energy = compute_xtb_energy(ase_atoms)
        all_energy[basename] = energy
        filetyper.append_json(all_energy, json_filename)
