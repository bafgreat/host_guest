#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
__credits__ = 'Dr. Matthew Addicoat'
###################################################################################################################
# A script to generate simple host-guest interaction complexes between a monomer and a host system.               #
# This script is a modification of the original Kick packages implemented by Dr. Matthew Addicoat.                #
# The original paper can be found in this paper 10.1002/jcc.21026                                                 #
###################################################################################################################

# import os
import operator
import random
from math import pi
from collections import namedtuple
# import glob
import numpy as np
# from ase.io import read, write
# from operator import attrgetter
from ase import geometry
from ase.data import *
from ase.neighborlist import *
from ase.atoms import Atoms
from ase.atom import Atom
from host_guest.energy import compute_sp
from host_guest.geometry.pore_analyser import radius_from_convexhall

# from mofbatteryml.io import filetyper


def check_coords(mol0, mol1):
    """
    Check the distance between two molecules
    """
    keep_coords = True
    tmp_mol = mol0 + mol1
    mol0_natoms = len(mol0)
    for a0 in range(mol0_natoms-1):
        for a1 in range(mol0_natoms, len(tmp_mol)):
            if tmp_mol[a0].symbol != "X" and tmp_mol[a1].symbol != "X":
                this_dist = tmp_mol.get_distance(a0, a1, mic=True)
                cutoff_dist = covalent_radii[tmp_mol[a0].number] + \
                    covalent_radii[tmp_mol[a1].number] + 0.4
                if this_dist < cutoff_dist:
                    keep_coords = False
                    return keep_coords

    com = mol1.get_center_of_mass()
    for axis in range(3):
        if mol0.pbc[axis]:
            if com[axis] > mol0.cell[axis, axis] or com[axis] < -mol0.cell[axis, axis]:
                keep_coords = False
                return keep_coords
    return keep_coords


def bond_radii():
    all_bond_radii = {symbol: [r, r-0.1, r-0.2]
                 for symbol, r in zip(chemical_symbols, covalent_radii)}
    mmbond = {'Cr': 4, 'Mo': 4, 'V': 3, 'Rh': 1, 'Ru': 2, 'W': 4,
              'Os': 3, 'Re': 4, 'Pt': 1, 'Tc': 4, 'Ir': 1, 'Pd': 1}
    return all_bond_radii, mmbond


def check_coords2(mol0, mol1):
    """
    Check distances between old and new molecules (Atoms objects)
    We'll add the molecules together such that we can use the mic facility.
    mol0 is the 'current good' mol and carries the pbc information, if any."""
    all_bond_radii, _ = bond_radii()
    keep_coords = False

    if mol0.get_pbc().any():
        distances = geometry.get_distances(mol0.positions, mol1.positions, np.array(
            mol0.cell.tolist()), pbc=mol0.get_pbc().any())

    elif mol1.get_pbc().any():
        distances = geometry.get_distances(mol0.positions, mol1.positions, np.array(
            mol1.cell.tolist()), pbc=mol1.get_pbc().any())
    else:
        distances = geometry.get_distances(mol0.positions, mol1.positions)

    shortest_distance = np.min(distances[1])
    i1, i2 = np.unravel_index(distances[1].argmin(), distances[1].shape)
    vdw_1 = all_bond_radii[mol0[i1].symbol][0]
    vdw_2 = all_bond_radii[mol1[i2].symbol][0]
    bond_cut_off = vdw_2 + vdw_1 + 1
    if shortest_distance > bond_cut_off:
        keep_coords = True
    return keep_coords


def molecule_diameter(new_atom):
    all_bond_radii, _ = bond_radii()
    dist_matrix = new_atom.get_all_distances(mic=True)
    final_matrix = np.triu(dist_matrix)
    i1, i2 = np.unravel_index(final_matrix.argmax(), final_matrix.shape)
    vdw_1 = all_bond_radii[new_atom[i1].symbol][0]
    vdw_2 = all_bond_radii[new_atom[i2].symbol][0]
    maxdim = final_matrix[i1, i2]+vdw_1 + vdw_2
    return maxdim


def Dock(host_system,  monomer, number_of_host=1, number_of_monomers=1, number_of_complexes=1, energy=True):
    """
    A fucntion to generate a host-guest system with a given number of monomers and number of hosts.
    parameter
    ----------
    host_system: ase.Atoms object
    monomer: ase.Atoms object
    number_of_host: int
    number_of_monomers: int
    number_of_complexes: int
    """
    print (f'energy in docker is{energy}')
    system_def = {}
    system_master = {}
    energy_dict = {}
    complex_molecules = {}
    k3_params = {'num_geoms': number_of_complexes, 'ncpus': 1, 'mem': '2GB', 'wall': 1,
                 'charge': 0, 'multiplicity': 1, 'program': 'dftb', 'order': True}
    if energy:
        energy_dict['monomer'] = compute_sp.compute_xtb_energy(monomer)
        energy_dict['host_system'] = compute_sp.compute_xtb_energy(host_system)

    SystemDict = namedtuple('SystemDict', 'number, frag_type')
    system_def['monomer'] = SystemDict(
        number=int(number_of_monomers), frag_type="R")

    system_master['monomer'] = monomer

    system_cell = host_system.get_cell()
    system_pbc = host_system.get_pbc()

    system_master['host_system'] = host_system
    framework = host_system
    if system_pbc.any():
        system_def['host_system'] = SystemDict(
        number=int(number_of_host), frag_type="P")
        k3_params["cell"] = system_cell
        k3_params["pbc"] = system_pbc
        k3_params["radius"] = molecule_diameter(framework)/2.0
    else:
        system_def['host_system'] = SystemDict(
        number=int(number_of_host), frag_type="L")
        k3_params["radius"] = radius_from_convexhall(framework.positions)
    # total_objects = sum([v.number for v in system_def.values()])

    frag_list = []
    tmp_sysdef = {}
    for k, v in system_def.items():
        print("H", k, v)
        tmp_sysdef[k] = int(v.number)
        if v.frag_type == "L" or v.frag_type == "P":
            frag_list.append(k)
            tmp_sysdef[k] -= 1

    if k3_params["order"]:
        while sum(tmp_sysdef.values()) != 0:
            for k, v in sorted(tmp_sysdef.items(), key=operator.itemgetter(1)):
                if v > 0:
                    frag_list.append(k)
                    tmp_sysdef[k] -= 1
    else:
        # for k, v in sorted(tmp_sysdef.iteritems(), key=operator.itemgetter(1)):
        #     for _ in range(v):
        #         frag_list.append(k)
        sorted_items = sorted(tmp_sysdef.items(), key=lambda x: x[1])
        for k, v in sorted_items:
            # Append the key to frag_list 'v' times
            frag_list.extend([k] * v)

    for ci in range(k3_params["num_geoms"]):
        # current = ci
        # this_file = "%04d" % current

        if "pbc" in k3_params:
            new_mol = Atoms(pbc=k3_params["pbc"], cell=k3_params["cell"])
            dummy_mol = Atoms(pbc=k3_params["pbc"], cell=k3_params["cell"])
        else:
            new_mol = Atoms()
            dummy_mol = Atoms()

        if k3_params["program"] == 'gaussian':
            new_mol.info["mem"] = k3_params["mem"]
            new_mol.info["ncpus"] = k3_params["ncpus"]
            # new_mol.info["filebase"] = filebase
            new_mol.info["charge"] = k3_params["charge"]
            new_mol.info["multiplicity"] = k3_params["multiplicity"]

        origin = Atom('X', position=(0, 0, 0))

        new_mol += origin
        dummy_mol += origin

        axes = [i for i in range(len(new_mol.get_pbc())) if new_mol.pbc[i]]
        new_mol.center(axis=axes)
        dummy_mol.center(axis=axes)
        new_mol.pop()

        for this_frag in frag_list:
            keep_coords = False
            attempts = 0
            max_attempts = 3000
            if system_def[this_frag].frag_type == "L" or system_def[this_frag].frag_type == "P":
                new_mol += system_master[this_frag]
                this_com = system_master[this_frag].get_center_of_mass()
                dummy_mol += Atom('X', position=this_com)
            else:
                this_origin = random.choice(dummy_mol).position
                while not keep_coords:
                    translation = np.random.random_sample(
                        3) * k3_params["radius"] * 2 - k3_params["radius"]
                    translation += this_origin
                    this_phi = random.uniform(-pi, pi)
                    this_theta = random.uniform(-pi, pi)
                    this_psi = random.uniform(-pi, pi)
                    tmp_mol = system_master[this_frag].copy()

                    if system_def[this_frag].frag_type != "A":
                        tmp_mol.euler_rotate(
                            center='COM', phi=this_phi*180/np.pi, theta=this_theta*180/np.pi, psi=this_psi*180/np.pi)

                    tmp_com = tmp_mol.get_center_of_mass()
                    tmp_mol.translate(translation - tmp_com)
                    attempts += 1
                    print ("attempts: ", attempts)
                    keep_coords = check_coords2(new_mol, tmp_mol)
                    if keep_coords:
                        new_mol += tmp_mol
                        tmp_com = tmp_mol.get_center_of_mass()
                        dummy_mol += Atom('X', position=tmp_com)
                        print ("complex build successfully after {} attempts.".format(attempts))
                        break
                    elif attempts == max_attempts:
                        print("Failed to add molecule after {} attempts. Moving to the next complex.".format(max_attempts))
                        new_mol = Atoms()
                        break
                        # this_origin = random.choice(dummy_mol).position
                        # attempts = 0
        if len(new_mol) > 0:
            host_guest_complex = 'complex_' + str(ci)
            if energy:
                new_mol_energy = compute_sp.compute_xtb_energy(new_mol)
                energy_dict[host_guest_complex] = new_mol_energy
                complex_molecules[host_guest_complex] = new_mol
                try:
                    free_energy = energy_dict[host_guest_complex]['energy_kcal_mol'] - (
                        number_of_monomers*energy_dict['monomer']['energy_kcal_mol'] + number_of_host*energy_dict['host_system']['energy_kcal_mol'])
                    energy_dict[host_guest_complex]['free_energy_kcal_mol'] = free_energy
                except Exception as e:
                    print(e)
            else:
                complex_molecules[host_guest_complex] = new_mol

    return energy_dict, complex_molecules
