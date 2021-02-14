"""
This is a suggestion for structuring your simulation code properly.
However, it is not set in stone. You may modify it if you feel like
you have a good reason to do so.
"""

import numpy as np

global positions_store, velocities_store

# initalizing self defined system parameters
num_atoms = 100    # amount of particles
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


# particle information
positions = init_position(N, M, dims)
velocities = init_velocity(N, temperature, dims)  # //np.zeros([N,dims])


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

    return


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
                rel_pos[i][j][k] = pos[j][k] - pos[i][k]

    rel_dist = np.zeros([len(pos), len(pos)])
    for i in range(0, len(rel_pos)):
        for j in range(0, len(rel_pos)):
            rel_dist[i][j] = np.math.sqrt(sum(i**2 for i in rel_pos[i][j]))

    return rel_pos, rel_dist


loc = init_position(num_atoms, box_dim, dim)
rel_loc = atomic_distances(loc, box_dim)
print(rel_loc[1][0])
print(len(rel_loc[1]))

# +
loc = init_position(num_atoms, box_dim, dim)
rel_loc = atomic_distances(loc, box_dim)
dUdt = np.zeros([num_atoms, num_atoms])

for i in range (0,len(rel_loc[1])): #particle i
    for j in range (0, len(rel_loc[1])): #particle i rel to j!=i
        if i != j:
            dUdt[i, j] = -24*EPSILON*((2*SIGMA**12/(rel_loc[1]**13)) - (SIGMA**6/rel_loc[1]**13))
        else:
            dUdt[i, j]

# -

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
    np.ndarray
        The net force acting on particle i due to all other particles
    """
    force = np.zeros(len(rel_dist))

    for i in range(0, len(rel_dist)):
        for j in range(0, len(rel_dist[i])):
            # does something with rel_dist
            if rel_dist[i][j] == 0.0:
                # do not include contributions to self
                continue

            # du = 4*EPSILON*((SIGMA/rel_dist[i][j])**12-(SIGMA/rel_dist[i][j])**6)
            du_dr = 4*EPSILON*(((SIGMA**12)*(-12)/(rel_dist[i][j]**13))-((SIGMA**6)*(-6)/(rel_dist[i][j])**7))
            force[i] -= du_dr

    return force


# x,y = atomic_distances(positions, M)
# print(lj_force(x,y))


def fcc_lattice(num_atoms, lat_const):
    """
    Initializes a system of atoms on an fcc lattice.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system
    lat_const : float
        The lattice constant for an fcc lattice

    Returns
    -------
    pos_vec : np.ndarray
        Array of particle coordinates
    """

    return


def kinetic_energy(vel):
    """
    Computes the kinetic energy of an atomic system.

    Parameters
    ----------
    vel: np.ndarray
        Velocity of particle

    Returns
    -------
    float
        The total kinetic energy of the system.
    """

    ke = 0;

    for i in range(0, len(vel)):
        ke += 0.5 * np.power(np.math.sqrt(sum(i ** 2 for i in vel[i])), 2.0)

    return ke


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
    num_atoms = len(rel_dist[0])
    # pot_e = np.zeros([num_atoms, num_atoms]) initializing causes problems, function still works
    for j in range(0, num_atoms):
        for i in range(0, num_atoms):
            if i != j:
                pot_e[i][j] = 4*EPSILON*((SIGMA/rel_dist[i][j])**12-(SIGMA/rel_dist[i][j])**6)
            else:
                pot_e[i][j] = 0
    pot_eparticle = np.sum(pot_e, axis=1)
    pot_total = np.sum(pot_eparticle)/2
    return pot_e, pot_eparticle, pot_total


def total_energy(vel, rel_dist):
    """
        Computes the total energy of an atomic system.

        Parameters
        ----------
        vel: np.ndarray
            Velocity of particle
        rel_dist : np.ndarray
            Relative particle distances as obtained from atomic_distances

        Returns
        -------
        float
            The total energy of the system.

        """

    return kinetic_energy(vel) + potential_energy(rel_dist)
