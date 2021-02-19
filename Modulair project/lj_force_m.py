from constants_m import *
import numpy as np

def lj_force(rel_pos, rel_dist):
    """
    Calculates the net forces on each atom.

    Parameters
    ----------
    rel_pos : np.ndarray
        Relative particle positions as obtained from atomic_distances
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    force : np.ndarray
        The force of atom j, on atom i. Where j are total-1 atoms.
    force_atom : np.ndarray
        The net force acting on particle i due to all other particles
        
    NOTE: THIS IS HOW INPUT CAN BE FOUND:
    loc = init_position(num_atoms, box_dim, dim)
    positions = atomic_distances(loc, box_dim)
    rel_dist = positions[1]
    rel_pos = positions[0]
    """
    dUdt = np.zeros([len(rel_dist), len(rel_dist)])
    force = np.zeros([len(rel_pos[1]), len(rel_pos[1]), len(rel_pos[0][0])])

    for i in range (0,len(rel_pos[1])): #particle i
        for j in range (0, len(rel_pos[1])): #particle i rel to j (!=i)
            if i != j:
                dUdt[i, j] = -24*EPSILON*((2*SIGMA**12/(rel_dist[i, j]**13)) - (SIGMA**6/rel_dist[i, j]**13))/(rel_dist[i, j]) 
            else:
                dUdt[i, j] = 0
    for i in range (0,len(rel_pos[1])): #particle i
        for j in range (0, len(rel_pos[1])): #particle i rel to j (!=i)
            force[i, j, :] = dUdt[i, j]*rel_pos[i, j, :]
    # while this looks horrible, and is horrible, it works. However, needs significant optimazation

    force_atom = np.sum(force, axis=1)

    return force, force_atom
