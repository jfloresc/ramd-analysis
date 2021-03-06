#!/usr/bin/env python
""" This script gets center of mass from a trajectory
by jose flores-canales 18.04.2018"""


from __future__ import print_function
import sys
import os
import optparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
#from sklearn.decomposition import PCA


def parse_cmdline(cmdlineargs):
    """parsing user defined parameters"""
    parser = optparse.OptionParser("Usage: python get_com.py --trajfile filename"
                                   "[options]")

    parser.add_option("-t", "--trajfile", action="store", dest="trajfile")
    
    opts, _ = parser.parse_args(cmdlineargs)
    traj_file = opts.trajfile
    if (traj_file is None): # or (topol_file is None) or (pdb_ref is None):
        parser.print_help()
        exit()
    return traj_file


def get_cog(traj):
    """get geometric center of the trajectory"""
    data = []
    n_frames = traj.n_frames 
    n_atoms = traj.n_atoms
    average = np.zeros((n_atoms, 3))
    for frame in xrange(n_frames):
        s_x, s_y, s_z = 0.0, 0.0, 0.0
        coord = traj.xyz[frame] 
        for atom_i in xrange(n_atoms):
            s_x += coord[atom_i, 0]
            s_y += coord[atom_i, 1]
            s_z += coord[atom_i, 2]
        s_x /= n_atoms
        s_y /= n_atoms
        s_z /= n_atoms
        data.append([s_x, s_y, s_z])
    return np.array(data)


def get_com_ligand(filename):
    """get center of mass from a trajectory"""

    print("Trajectory: ", filename)
    topol_name = 'topol.prmtop'
    pdb_ref = 'start_conf.pdb'
    T1 = md.load(filename, top=topol_name)
    topol = md.load(pdb_ref)
    topology = topol.topology
    selection = 'resname GLA'
    atom_to_keep = topology.select(selection)
    print("Number of frames in trajectory:", T1.n_frames)
    T1.restrict_atoms(atom_to_keep)
    print("Number of atoms in selection:", T1.n_atoms)
    coord = md.compute_center_of_mass(T1)
    coord_cog = get_cog(T1)
    print("Dimensions of COM coordinates:", coord.shape, "\n")

    basename = os.path.basename(filename)
    out = basename.split('.')
    output_file = out[0] + '_COM.dat'
    np.savetxt(output_file, coord, fmt='%3.4f', delimiter='\t') # Save array	
    output_file = out[0] + '_COG.dat'
    np.savetxt(output_file, coord_cog, fmt='%3.4f', delimiter='\t') # Save array	


###########################################################
# __main__
###########################################################
if __name__ == "__main__":
    TRAJ_FILENAME = parse_cmdline(sys.argv[1:]) 
    get_com_ligand(TRAJ_FILENAME)
