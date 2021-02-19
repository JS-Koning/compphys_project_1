from constants_m import *
import numpy as np


def kinetic_energy(vel):
    """
    Computes the kinetic energy of an atomic system.

    Parameters
    ----------
    vel: np.ndarray
        Velocity of particle

    Returns
    -------
    ke : float
        The total kinetic energy of the system.
    """

    ke = 0;

    for i in range(0, len(vel)):
        ke += 0.5 * ARG_MASS * np.power(np.math.sqrt(sum(i ** 2 for i in vel[i])), 2.0)

    if dimless:
        ke *= dimless_energy

    return ke
