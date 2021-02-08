"""
auto_cvs.py
Implement automatic CV definitions by building every 2nd, 3rd, and 4th order CV (bonds, angles, dihedrals) along bonded
atoms in the given topology file.
"""

import sys
import os
import mdtraj
import copy
import numpy
import collections
import itertools

import argparse

def main(settings):
    """
    Build and return a "cvs" list based on the topology file and other attributes of settings.

    Builds a CV definition for every 2nd, 3rd, and 4th order CV (bonds, angles, and dihedrals) along bonded atoms within
    settings.auto_cvs_radius (in angstroms) of each atom in either commitment basin definition in the given topology
    file, returning a cvs list of the proper format for settings.cvs. Distances are based on the first coordinate file
    in settings.initial_coordinates.

    Also builds a text file in settings.working_directory describing each CV in the same order as they appear in cvs.

    Parameters
    ----------
    settings : argparse.Namespace
        Global settings object, including minimum attributes "working_directory", "topology", "auto_cvs_radius",
        "commit_fwd", "commit_bwd", and "initial_coordinates"

    Returns
    -------
    cvs : list
        List of CV definition strings

    """

    if settings.auto_cvs_type not in ['pytraj', 'mdtraj']:
        raise RuntimeError('unsupported auto_cvs_type setting: ' + settings.auto_cvs_type)

    # Initialize cvs list to eventually return, plus descriptions to pass to the user
    cvs = []
    descriptions = []

    # Load topology and coordinates as mdtraj "trajectory" object
    if settings.initial_coordinates[0] == '':
        raise RuntimeError('it appears that you forgot to provide an initial coordinate file for use with auto_cvs.py.')
    mtraj = mdtraj.load(settings.initial_coordinates[0], top=settings.topology)

    # Identify atoms within settings.auto_cvs_radius of each atom in either commitment basin definition
    commit_atoms = []   # atom indices involved in commitment basin definitions, as integers
    commit_atoms += settings.commit_fwd[0]
    commit_atoms += settings.commit_fwd[1]
    commit_atoms += settings.commit_bwd[0]
    commit_atoms += settings.commit_bwd[1]
    commit_atoms = list(set(commit_atoms))    # remove duplicates
    commit_atoms = [item - 1 for item in commit_atoms]      # fix off-by-one error

    # Compute neighbors with handy mdtraj function; divide radius by 10 to convert nm to Å
    neighbors = list(mdtraj.compute_neighbors(mtraj, settings.auto_cvs_radius / 10, query_indices=commit_atoms)[0])
    neighbors += commit_atoms   # include commit_atoms in neighbors for our purposes

    # Remove waters from neighbors, if necessary
    if settings.auto_cvs_exclude_water:
        temp = copy.deepcopy(neighbors)
        waters = mtraj.topology.select('water')
        for atom_index in neighbors:
            if atom_index in waters:
                temp.remove(atom_index)
        neighbors = copy.deepcopy(temp)

    # Assemble list of each 2nd order term
    bonds = []
    table, all_bonds = mtraj.topology.to_dataframe()
    formatted_all_bonds = [[int(item[0]), int(item[1])] for item in all_bonds]
    for atom_index in neighbors:
        bonds += [item for item in formatted_all_bonds if atom_index in item]
    unq_lst = collections.OrderedDict()
    for item in bonds:
        unq_lst.setdefault(frozenset(item), []).append(item)
    bonds = list(map(list, unq_lst.keys()))

    # Add every pair of atoms in commit_fwd and/or commit_bwd to the list of "bonds" if not already present
    for first_index in commit_atoms:
        for second_index in [item for item in commit_atoms if not item == first_index]:
            if not [first_index, second_index] in bonds and not [second_index, first_index] in bonds:
                bonds.append([first_index, second_index])

    # Append code to obtain each bond length as a string to cvs
    for bond in bonds:
        if settings.auto_cvs_type == 'mdtraj':
            cv_str = 'mdtraj.compute_distances(mtraj, numpy.array([' + str(bond) + ']))[0][0] * 10'
        else:
            cv_str = 'pytraj.distance(traj, \'@' + ' @'.join([str(int(atom + 1)) for atom in bond]) + '\')[0]'
        cvs.append(cv_str)
        descriptions.append('distance between atoms ' + str([str(int(atom + 1)) for atom in bond]))

    # Repeat previous process for angles; only difference is added logic to find bonded triplets since there is no
    # equivalent to topology._bonds for angles
    angles = []
    for first_index in neighbors:
        for second_index in [item for item in neighbors if not item == first_index]:
            for third_index in [item for item in neighbors if not item == first_index and not item == second_index]:
                if {first_index, third_index} in [set(item) for item in bonds] and {second_index, third_index} in [set(item) for item in bonds]:
                    if not [first_index, third_index, second_index] in angles and not [second_index, third_index, first_index] in angles:  # prevent duplicates
                        angles.append([first_index, third_index, second_index])

    # Append code to obtain each angle as a string to cvs
    for angle in angles:
        if settings.auto_cvs_type == 'mdtraj':
            cv_str = 'mdtraj.compute_angles(mtraj, numpy.array([' + str(angle) + ']))[0][0] * 180 / numpy.pi'
        else:
            cv_str = 'pytraj.angle(traj, \'@' + ' @'.join([str(int(atom + 1)) for atom in angle]) + '\')[0]'
        cvs.append(cv_str)
        descriptions.append('angle between atoms ' + str([str(int(atom + 1)) for atom in angle]))

    # Repeat for dihedrals
    dihedrals = []
    for angle in angles:
        for bond in bonds:
            if angle[0] in bond and not angle[1] in bond:
                temp = copy.deepcopy(bond)
                temp.remove(angle[0])
                if not reversed(temp + angle) in dihedrals and not temp in angle:
                    dihedrals.append(temp + angle)
            if angle[2] in bond and not angle[1] in bond:
                temp = copy.deepcopy(bond)
                temp.remove(angle[2])
                if not reversed(angle + temp) in dihedrals and not temp in angle:
                    dihedrals.append(angle + temp)

    # Append code to obtain each dihedral as a string to cvs
    for dihedral in dihedrals:
        if settings.auto_cvs_type == 'mdtraj':
            cv_str = 'mdtraj.compute_dihedrals(mtraj, numpy.array([' + str(dihedral) + ']))[0][0] * 180 / numpy.pi'
        else:
            cv_str = 'pytraj.dihedral(traj, \'@' + ' @'.join([str(int(atom + 1)) for atom in dihedral]) + '\')[0]'
        cvs.append(cv_str)
        descriptions.append('dihedral between atoms ' + str([str(int(atom + 1)) for atom in dihedral]))

    # And finally, difference-of-distance terms for atoms that appear in bonds with at least two other atoms in the
    # commitment definitions
    partners = []
    for atom_index in commit_atoms:
        this_partners = []
        for local_index in range(len(settings.commit_fwd[0])):
            if settings.commit_fwd[0][local_index] - 1 == atom_index and not settings.commit_fwd[1][local_index] - 1 in this_partners:
                this_partners += [settings.commit_fwd[1][local_index] - 1]
        for local_index in range(len(settings.commit_fwd[1])):
            if settings.commit_fwd[1][local_index] - 1 == atom_index and not settings.commit_fwd[0][local_index] - 1 in this_partners:
                this_partners += [settings.commit_fwd[0][local_index] - 1]
        for local_index in range(len(settings.commit_bwd[0])):
            if settings.commit_bwd[0][local_index] - 1 == atom_index and not settings.commit_bwd[1][local_index] - 1 in this_partners:
                this_partners += [settings.commit_bwd[1][local_index] - 1]
        for local_index in range(len(settings.commit_bwd[1])):
            if settings.commit_bwd[1][local_index] - 1 == atom_index and not settings.commit_bwd[0][local_index] - 1 in this_partners:
                this_partners += [settings.commit_bwd[0][local_index] - 1]
        partners.append(this_partners)
    for partners_index in range(len(partners)):
        if len(partners[partners_index]) < 1:   # just a quick fail-safe
            raise RuntimeError('Internal error in auto_cvs.py: commit_atom ' + str(commit_atoms[partners_index]) + ' is'
                               ' without a partner.')
        if len(partners[partners_index]) > 1:   # add difference-of-distance terms for this atom
            for combination in itertools.combinations(partners[partners_index], 2):
                if settings.auto_cvs_type == 'mdtraj':
                    cv_str = '(mdtraj.compute_distances(mtraj, numpy.array([[' + str(commit_atoms[partners_index]) + ', ' + str(combination[0]) + ']]))[0][0] * 10) - ' +\
                             '(mdtraj.compute_distances(mtraj, numpy.array([[' + str(commit_atoms[partners_index]) + ', ' + str(combination[1]) + ']]))[0][0] * 10)'
                else:
                    cv_str = 'pytraj.distance(traj, \'@' + str(int(commit_atoms[partners_index] + 1)) + ' @' + str(int(combination[0] + 1)) + '\')[0] - ' + \
                             'pytraj.distance(traj, \'@' + str(int(commit_atoms[partners_index] + 1)) + ' @' + str(int(combination[1] + 1)) + '\')[0]'
                cvs.append(cv_str)
                descriptions.append('difference of distances between atoms [' + str(int(commit_atoms[partners_index] + 1)) + ', ' + str(int(combination[0] + 1)) + '] and [' + str(int(commit_atoms[partners_index] + 1)) + ', ' + str(int(combination[1] + 1)) + ']')

    # Now just create the output text document and return
    if not os.path.exists(settings.working_directory):
        os.mkdir(settings.working_directory)
    open(settings.working_directory + '/cvs.txt', 'w').write('CV name: description of CV (atoms are indexed from 1); code to evaluate CV\n')
    cv_index = 0    # initialize CV index
    with open(settings.working_directory + '/cvs.txt', 'a') as f:
        for cv in cvs:
            cv_index += 1
            f.write('CV' + str(cv_index) + ': ' + descriptions[cv_index - 1] + '; ' + cv + '\n')
        if settings.cvs:
            for cv in settings.cvs:
                if not cv == '':    # this is, if this CV is actually defined
                    cv_index += 1
                    cvs.append(cv)
                    f.write('CV' + str(cv_index) + ': user-defined CV; ' + cv + '\n')

    return cvs
