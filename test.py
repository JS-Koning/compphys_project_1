import numpy as np
import matplotlib.pyplot as plt
global positions_store, velocities_store
# initalizing self defined system parameters
num_atoms = 9  # amount of particles
dim = 3  # dimensions
box_dim = 10  # meters; bounding box dimension
dt = 1e-4  # s; stepsize
steps = 10000  # amount of steps
dimless = True  # use dimensionless units
periodic = True  # use periodicity
verlet = True  # use Verlet's algorithm
# Parameters physical, supplied by course, or related to Argon
temp = 119.8  # Kelvin
KB = 1.38064852e-23  # m^2*kg/s^2/K
SIGMA = 3.405e-10  # meter
EPSILON = temp * KB  # depth of potential well/dispersion energy
N_b = 6.02214076e23  # Avagadros number; 1/mol
R = 8.31446261815324  # J/K/mole; universal gas constant
ARG_UMASS = 39.95  # u; atomic mass of argon
ARG_MMASS = ARG_UMASS / 1000  # kg/mol; mole mass of argon
ARG_MASS = ARG_UMASS * 1.6605e-27  # Kg mass of a single atom in Kg
# conversion values for dimensionless units
dimless_time = 1.0 / np.math.sqrt((ARG_MASS * SIGMA ** 2 / EPSILON))  # s; dimensionless time
dimless_energy = 1.0 / EPSILON  # J; dimensionless energy
dimless_distance = 1.0 / SIGMA  # m; dimensionless distance
dimless_velocity = 1.0 / np.math.sqrt(EPSILON / (ARG_MASS))  # m/s; dimensionless velocity


def init_velocity(num_atoms, temp, dim):
    """
    Initializes the system with Gaussian distributed velocities. This
    init_velocity is loosely based on 3D system, however it will output 2D just
    fine, although more pertubated. Note, relativity is ignored.
    For more information, please visit articles related to Maxwell Boltzmann
    distributions.
    https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
    
    Parameters
    ----------
    num_atoms : int
        The number of particles in the system.
    temp : float
        The (unitless) temperature of the system.
    dim : int
        The dimensions of the system.

    Returns
    -------
    vel_vec : np.ndarray
        Array of particle velocities
    """
    # define most probable speed (vel_p), then use this to find mean speed, using Maxwell-Boltzmann distribution
    vel_p = np.sqrt(2 * KB * temp / ARG_MASS)

    if dimless:
        vel_p *= dimless_velocity

    vel_mean = 2 * vel_p / np.sqrt(np.pi)
    # define the standard deviation, assuming the standard deviation: std = sqrt(meansquarevelocity^2 - mean velocity^2) 
    # again, using Maxwell-Boltzmann distributions
    vel_msq = (3 * vel_p ** 2) / 2
    vel_std = vel_msq - (vel_mean ** 2)
    # find the final distribution
    vel_vec = np.random.normal(vel_mean, vel_std, (num_atoms, dim))
    return vel_vec


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
        if dimless:
            velsquared = np.power(vel, 2.0)
            vel_summed = np.sum(velsquared, axis=1)
            vel_abs = np.power(vel_summed, 0.5) #the total velocity, of 1 particle stored in an array for each particle. Since a bug was present, This is rewritten in, over simplified steps.
            ke_part = 0.5 * np.power(vel_abs, 2)
            ke_total = np.sum(ke_part)

        else:
            ke += 0.5 * ARG_MASS * np.power(np.math.sqrt(sum(i ** 2 for i in vel[i])), 2.0)

    return ke_part, ke_total


b = np.array([3,0,0])
print(b)
#a = kinetic_energy(b)

a = init_velocity(3,273,3)
b = a + 1
c = np.append(a, b, axis=0)
print(c)

kinetic_energy(c)




