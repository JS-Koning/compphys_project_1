"""
This is a suggestion for structuring your simulation code properly.
However, it is not set in stone. You may modify it if you feel like
you have a good reason to do so.
"""

import numpy as np

global positions_store, velocities_store

# initalizing self defined system parameters
num_atoms = 50    # amount of particles
dim = 2            # dimensions
box_dim = 10      # meters; bounding box dimension
dt = 0.01          # s; stepsize
steps = 100        # amount of steps

# Parameters physical, supplied by course, or related to Argon
temp = 119.8                # Kelvin
KB = 1.38064852e-23         # m^2*kg/s^2/K
SIGMA = 3.405e-10           # meter
EPSILON = temp * KB         # depth of potential well/dispersion energy
N_b = 6.02214076e23         # Avagadros number
R = 8.31446261815324        # J/K/mole; universal gas constant
ARG_UMASS = 39.95           # u; atomic mass of argon
ARG_MMASS = ARG_UMASS/1000  # kg/mol; mole mass of argon



def init_velocity(num_atoms, temp, dim):
    
    
    """
    Initializes the system with Gaussian distributed velocities.

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
    vel_p = np.sqrt(2*R*temp/ARG_MMASS)
    vel_mean = 2*vel_p/np.sqrt(np.pi)
    # define the standard deviation, assuming the standard deviation: std = sqrt(meansquarevelocity^2 - mean velocity^2) 
    # again, using Maxwell-Boltzmann distributions
    vel_msq = (3*vel_p**2)/2
    vel_std = vel_msq-(vel_mean**2)
    # find the final distribution
    vel_vec = np.random.normal(vel_mean, vel_std, (num_atoms, dim))
    return vel_vec


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
    pos_vec = np.random.random((num_atoms, dim)) * box_dim
    return pos_vec

def simulate(init_pos, init_vel, num_tsteps, timestep, box_dim):
    """
    Molecular dynamics simulation using the Euler or Verlet's algorithms
    to integrate the equations of motion. Calculates energies and other
    observables at each timestep.

    Parameters
    ----------
    init_pos : np.ndarray
        The initial positions of the atoms in Cartesian space
    init_vel : np.ndarray
        The initial velocities of the atoms in Cartesian space
    num_tsteps : int
        The total number of simulation steps
    timestep : float
        Duration of a single simulation step
    box_dim : float
        Dimensions of the simulation box

    Returns
    -------
    Any quantities or observables that you wish to study.
    """
     # first initialize matrix starting with initial velocity and position
    pos_steps = np.zeros((num_tsteps, num_atoms, dim))
    vel_steps = np.zeros((num_tsteps, num_atoms, dim))
    init_pos = init_position(num_atoms, box_dim, dim)
    init_vel = init_velocity(num_atoms, box_dim, dim)
    # everything need fixing <start>
    init_rel_pos = atomic_distances(pos, box_dim, dim)
    init_force = lj_force(rel_pos, rel_dist, dim)
    # everything need fixing </end>
    pos_steps[0, :, :] = init_pos
    vel_steps[0, :, :] = init_vel
    for i in range steps:
        pos_steps[i+1, :, :] = pos_steps[i, :, :] + vel_steps[i, :, :]*timestep
        vel_steps[i+1, :, :] = vel_steps[i, :, :] + force[i, :, :]*timestep/ARG_UMASS # NOTE, this mass is not yet correct! 
        return pos_steps, vel_steps
    return pos_steps, vel_steps


loc = init_position(num_atoms, box_dim, dim)

print(loc)

print(len(loc))
print(loc.shape)
qq = np.zeros(((num_atoms, dim), steps))
print(qq.shape)
print(qq)


y = np.arange(3, 6)
print(y)

loc_steps = np.zeros((dim, len(loc), steps))
print(loc_steps)

loc_steps.shape

qq = np.zeros((50,2))
qq.shape
print(qq)


loc_steps[:,:,0] = loc[:,:]
print(loc_steps)
loc_steps.shape
