import numpy as np
from constants_m import *


def atomic_distances(pos, box_dim):
    """
    Calculates relative positions and distances between particles.

    parameters
    ----------
    pos : np.ndarray
        The positions of the particles in cartesian space
    box_dim : float
        The dimension of the simulation box

    returns
    -------
    rel_pos : np.ndarray
        Relative positions of particles
    rel_dist : np.ndarray
        The distance between particles
    """

    dim = len(pos[0])
    # NOTE: includes rel_dist/rel_pos to itself (= 0 / [0.0, 0.0])
    rel_pos = np.zeros([len(pos), len(pos), dim])

    for i in range(0, len(pos)):
        for j in range(0, len(pos)):
            for k in range(0, dim):
                dis = pos[j][k] - pos[i][k]
                if periodic:
                    if dimless:
                        if dis > (box_dim * dimless_distance * 0.5):
                            dis = dis - box_dim * dimless_distance
                        if dis <= -(box_dim * dimless_distance * 0.5):
                            dis = dis + box_dim * dimless_distance
                    else:
                        if dis > (box_dim * 0.5):
                            dis = dis - box_dim
                        if dis <= -(box_dim * 0.5):
                            dis = dis + box_dim
                rel_pos[i][j][k] = dis

    rel_dist = np.zeros([len(pos), len(pos)])
    for i in range(0, len(rel_pos)):
        for j in range(0, len(rel_pos)):
            rel_dist[i][j] = np.math.sqrt(sum(i**2 for i in rel_pos[i][j]))

    return rel_pos, rel_dist

