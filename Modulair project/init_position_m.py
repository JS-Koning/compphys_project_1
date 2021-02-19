from constants_m import *
import numpy as np


def init_position(num_atoms, box_dim, dim):
    """
    Initializes the system with random positions.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system.
    box_dim : float
        The dimension of the simulation box
    dim : int
        The dimensions of the system.

    Returns
    -------
    pos_vec : np.ndarray
        Array of particle positions
    """

    random = np.random.random((num_atoms, dim))

    # to test close to boundary
    #random *= 0.001

    pos_vec = random * box_dim

    if dimless:
        pos_vec *= dimless_distance

    return pos_vec
