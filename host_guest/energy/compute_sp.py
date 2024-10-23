import os
import subprocess
import shutil
import random
from host_guest.io import coords_library, filetyper

class ComputeEnergy():
    """
    A class to compute energy for a given system. There are two methods to
    compute the energies. This includes 'xtb' and 'machine learning potential
    from orb-models.

    Method 1: xtb
    ---------------

    Method 2: Machine Learning Potential from Orb-Model
    -----------------------------------------------------
    The



    """
    def __init__(self, ase_atoms, method='xtb'):
        self.ase_atoms = ase_atoms
        self.method = method



    def use_xtb(self):
        """
        A function to determine which method to use. The function
        checks whether xtb is installed. If method is set to xtb and
        xtb is installed then the model will use xtb. Otherwise the
        model will use the default **orb_d3_v1**. Other methods that
        be set are
        - ***xtb**
        - **orb_v1**
        - ***xtb**
        - **orb_mptraj_only_v1**
        - **orb_d3_{sm,xs}_v1**
        """
        xtb_method = False
        try:
            result = subprocess.run(["xtb", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0 and self.method=='xtb':
                xtb_method = True
        except FileNotFoundError:
            print ('xtb is not installed please run ')
            pass
        return xtb_method


    def compute_xtb_energy(self):
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
        coords_library.write_ase_atoms(self.ase_atoms, tmp_in)
        os.system(f'xtb --sp --gfn 2 --tblite --spinpol {tmp_in} > {tmp_out}')
        energy = read_xtb_energy(tmp_out)
        os.chdir(base_dir)
        if os.path.exists(result_folder):
            shutil.rmtree(result_folder)
        return energy


    def read_xtb_energy(self, filename):
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


    def compute_energy_from_orb_model(self):
        """
        A function that computes the energy of a given ase atom using
        a machine learning potential from orb-models.
        """
        # orbff = pretrained.orb_d3_v1()
        ase_graph = coords_library.ase_graph(self.ase_atoms)
        # result = orbff.predict(ase_graph)
        # print (result)




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
