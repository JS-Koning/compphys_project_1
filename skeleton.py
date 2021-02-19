"""
This is a suggestion for structuring your simulation code properly.
However, it is not set in stone. You may modify it if you feel like
you have a good reason to do so.
"""

import numpy as np
import matplotlib.pyplot as plt

global positions_store, velocities_store

# initalizing self defined system parameters
num_atoms = 2               # amount of particles
dim = 2                     # dimensions
box_dim = 10                # meters; bounding box dimension
dt = 0.01                   # s; stepsize
steps = 100                 # amount of steps
dimless = True              # use dimensionless units
periodic = True             # use periodicity

# Parameters physical, supplied by course, or related to Argon
temp = 119.8                # Kelvin
KB = 1.38064852e-23         # m^2*kg/s^2/K
SIGMA = 3.405e-10           # meter
EPSILON = temp * KB         # depth of potential well/dispersion energy
N_b = 6.02214076e23         # Avagadros number; 1/mol
R = 8.31446261815324        # J/K/mole; universal gas constant
ARG_UMASS = 39.95           # u; atomic mass of argon
ARG_MMASS = ARG_UMASS/1000  # kg/mol; mole mass of argon
ARG_MASS = ARG_UMASS*1.6605e-27 #Kg mass of a single atom in Kg

# conversion values for dimensionless units
dimless_time = 1.0 / np.math.sqrt((ARG_MASS*SIGMA**2/EPSILON)) # s; dimensionless time
dimless_energy = 1.0 / EPSILON                                      # J; dimensionless energy
dimless_distance = 1.0 / SIGMA                                      # m; dimensionless distance
dimless_velocity = 1.0 / np.math.sqrt(EPSILON/(ARG_MASS))      # m/s; dimensionless velocity

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
    vel_p = np.sqrt(2*KB*temp/ARG_MASS)

    if dimless:
        vel_p *= dimless_velocity

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

    random = np.random.random((num_atoms, dim))

    # to test close to boundary
    #random *= 0.001

    pos_vec = random * box_dim

    if dimless:
        pos_vec *= dimless_distance

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

    pos_steps[0, :, :] = init_pos
    vel_steps[0, :, :] = init_vel
    for i in range(num_tsteps-1):
        pos = pos_steps[i, :, :]

        # make sure it's inside box dimension -> modulus gives periodicity
        if periodic:
            if dimless:
                pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep*dimless_time) % box_dim*dimless_distance
            else:
                pos_steps[i+1, :, :] = (pos + vel_steps[i, :, :] * timestep ) % box_dim
        else:
            if dimless:
                pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep*dimless_time)
            else:
                pos_steps[i+1, :, :] = (pos + vel_steps[i, :, :] * timestep )

        rel_pos = atomic_distances(pos, box_dim)[0]
        rel_dis = atomic_distances(pos, box_dim)[1]
        force = lj_force(rel_pos, rel_dis)[1]

        if dimless:
            vel_steps[i + 1, :, :] = vel_steps[i, :, :] + force * timestep*dimless_time / (ARG_MMASS / N_b)
        else:
            vel_steps[i+1, :, :] = vel_steps[i, :, :] + force * timestep / (ARG_MMASS/N_b)

    global positions_store
    positions_store = pos_steps
    global velocities_store
    velocities_store = vel_steps

    return pos_steps, vel_steps


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
    ke : float
        The total kinetic energy of the system.
    """

    ke = 0;

    for i in range(0, len(vel)):
        ke += 0.5 * (ARG_MMASS/N_b) * np.power(np.math.sqrt(sum(i ** 2 for i in vel[i])), 2.0)

    if dimless:
        ke *= dimless_energy

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
        -------
        This is simply potential_energy[2]+kinetic_energy[0]

        """

    return kinetic_energy(vel) + potential_energy(rel_dist)[2]


def main():
    """"
        Beginning of program
    """
    init_pos = init_position(num_atoms, box_dim, dim)
    init_vel = init_velocity(num_atoms, box_dim, dim)

    simulate(init_pos, init_vel, steps, dt, box_dim)

    print("Test if the total energy is conserved")
    pos1 = positions_store[0, :, :]
    pos2 = positions_store[steps-1, :, :]

    vel1 = velocities_store[0, :, :]
    vel2 = velocities_store[steps-1, :, :]

    r_pos1 = atomic_distances(pos1, box_dim)
    r_pos2 = atomic_distances(pos2, box_dim)

    print("Initial total energy: " + str(total_energy(vel1, r_pos1[1])))
    print("Final total energy:   " + str(total_energy(vel2, r_pos2[1])))
    print("Delta total energy:   " + str(total_energy(vel2, r_pos2[1])-total_energy(vel1, r_pos1[1])))

    if num_atoms == 2:
        print("Plot inter-atom distance over time")
        if dimless:
            distances = [np.max(atomic_distances(positions_store[x, :, :],box_dim)[1])/dimless_distance
                         for x in range(steps)]
        else:
            distances = [np.max(atomic_distances(positions_store[x, :, :], box_dim)[1]) for x in range(steps)]
        times = np.linspace(0, dt*steps, steps)
        plt.plot(times, distances)
        plt.xlabel('Time (s)')
        plt.ylabel('Distance (m)')
        plt.show()

        print("Print energy levels over time")
        if dimless:
            energies = [(kinetic_energy(velocities_store[x, :, :])/dimless_energy,
                         potential_energy(atomic_distances(positions_store[x, :, :],box_dim)[1])[2]/dimless_energy,
                        total_energy(velocities_store[x, :, :],atomic_distances(positions_store[x, :, :],box_dim)[1])/dimless_energy)
                        for x in range(steps)]
        else:
            energies = [kinetic_energy(velocities_store[x, :, :]) for x in range(steps)]
        # times = np.linspace(0, dt*steps, steps)
        plt.plot(times, energies)
        plt.xlabel('Time (s)')
        plt.ylabel('Energy (J)')
        plt.show()


main()

