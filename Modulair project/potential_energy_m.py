import numpy as np
from constants_m import *

def potential_energy(rel_dist):
    """
    Computes the potential energy of an atomic system.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances
    NOTE!
    pos = init_position(num_atoms, box_dim, dim)
    rel_dist = atomic_distances(pos, box_dim)[1]
    !
    Returns
    -------
    pot_e : float
        The potential energy of a single atom, of each other atom.
    pot_etotal : float
        The potential energy of the atom of all other atoms
        NOTE: RETRIEVE BY print(pot[1])
    pot_total : float
        The total potential energy of the system
    """
    num_atoms1 = len(rel_dist[0])
    pot_e = np.zeros([num_atoms1, num_atoms1])
    for j in range(0, num_atoms1):
        for i in range(0, num_atoms1):
            if i != j:
                pot_e[i][j] = 4*EPSILON*((SIGMA/rel_dist[i][j])**12-(SIGMA/rel_dist[i][j])**6)
            else:
                pot_e[i][j] = 0
    pot_e_particle = np.sum(pot_e, axis=1)
    pot_total = np.sum(pot_e_particle)/2

    if dimless:
        pot_e *= dimless_energy
        pot_e_particle *= dimless_energy
        pot_total *= dimless_energy

    return pot_e, pot_e_particle, pot_total