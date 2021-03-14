# -*- coding: utf-8 -*-
"""
This is a suggestion for structuring your simulation code properly.
However, it is not set in stone. You may modify it if you feel like
you have a good reason to do so.
"""

import numpy as np
import random
import matplotlib.pyplot as plt

global positions_store, velocities_store

# initalizing self defined system parameters
num_atoms = 4  # amount of particles
dim = 3  # dimensions
box_dim = 1.547 #10  # meters; bounding box dimension
dt = 1e-4  # s; stepsize
steps = 30000  # amount of steps
dimless = True  # use dimensionless units
periodic = True  # use periodicity
verlet = True  # use Verlet's algorithm
rescaling = True # use Temperature rescaling
rescaling_mode = 1 # 0 = kin-NRG-based | temp-based
rescaling_delta = 0.0027 # delta for activation of rescaling
rescaling_timesteps = steps / 30 # timesteps interval for rescaling check
rescaling_max_timesteps = steps/2 # max timesteps for rescaling

# Parameters physical, supplied by course, or related to Argon
temp = 70  # Kelvin
TEMP = 119.8  # Kelvin
KB = 1.38064852e-23  # m^2*kg/s^2/K
SIGMA = 3.405e-10  # meter
EPSILON = TEMP * KB  # depth of potential well/dispersion energy
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


def init_velocity_old(num_atoms, temp, dim):
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
    # define most probable speed (vel_p)
    vel_p = np.sqrt(2 * KB * temp / ARG_MASS)

    if dimless:
        vel_p *= dimless_velocity

    # then use this to find mean speed, using Maxwell-Boltzmann distribution
    vel_mean = 2 * vel_p / np.sqrt(np.pi)

    # define the standard deviation, assuming the standard deviation: std = sqrt(meansquarevelocity^2 - mean velocity^2)
    # again, using Maxwell-Boltzmann distributions
    vel_msq = (3 * vel_p ** 2) / 2
    vel_std = vel_msq - (vel_mean ** 2)

    # find the distribution
    vel_vec = np.random.normal(vel_mean, vel_std, (num_atoms, dim))

    # get the magnitudes of the velocities
    vel_mag = [np.sqrt(np.sum([v[i]**2 for i in range(dim)])) for v in vel_vec]

    # rescale the magnitudes to match the vel_mean speed
    vel_vec *= vel_mean / np.mean(vel_mag)

    # create random negativity
    for v in range(num_atoms):
        for i in range(dim):
            # either *1 or *-1
            vel_vec[v,i] *= (1 - 2*random.getrandbits(1))

    # remove the mean velocity to keep 0 mean (no drift velocity)
    vel_vec -= np.mean(vel_vec)

    return vel_vec


def init_position(num_atoms, box_dim, dim):
    """
    Initializes the system with random positions.
    This does not require non dimensionalization scaling, since it is not
    based on physical parameters.

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
    # random *= 0.001

    pos_vec = random * box_dim

    return pos_vec


def simulate_old(init_pos, init_vel, num_tsteps, timestep, box_dim):
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
    for i in range(num_tsteps - 1):
        pos = pos_steps[i, :, :]

        # make sure it's inside box dimension -> modulus gives periodicity
        if periodic:
            if dimless:
                pos_steps[i + 1, :, :] = (pos + vel_steps[i, :,
                                                :] * timestep * dimless_time) % box_dim * dimless_distance
            else:
                pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep) % box_dim
        else:
            if dimless:
                pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep * dimless_time)
            else:
                pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep)

        rel_pos = atomic_distances(pos, box_dim)[0]
        rel_dis = atomic_distances(pos, box_dim)[1]
        force = lj_force(rel_pos, rel_dis)[1]

        if dimless:
            vel_steps[i + 1, :, :] = vel_steps[i, :, :] + force * timestep * dimless_time / ARG_MASS
        else:
            vel_steps[i + 1, :, :] = vel_steps[i, :, :] + force * timestep / ARG_MASS

    global positions_store
    positions_store = pos_steps
    global velocities_store
    velocities_store = vel_steps

    return pos_steps, vel_steps


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
    if verlet:
        print("Simulating using Verlet's algorithm")
    else:
        print("Simulating using Euler's algorithm")

    # first initialize matrix starting with initial velocity and position
    pos_steps = np.zeros((num_tsteps, num_atoms, dim))
    vel_steps = np.zeros((num_tsteps, num_atoms, dim))

    pos_steps[0, :, :] = init_pos
    vel_steps[0, :, :] = init_vel

    rescale_counter = 0
    rescale_max = 1.0
    rescale_min = 1.0

    for i in range(num_tsteps - 1):
        pos = pos_steps[i, :, :]

        if verlet:
            # Verlet velocity algorithm
            rel_pos, rel_dis = atomic_distances(pos, box_dim)
            force = lj_force(rel_pos, rel_dis)[1]
            #print(force)
            if periodic:
                # make sure it's inside box dimension -> modulus gives periodicity
                if dimless:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep + (
                            timestep ** 2) * force / 2) % box_dim
                else:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep + (
                            timestep ** 2) * force / 2) % box_dim
            else:
                if dimless:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep + (timestep ** 2) * force / 2)
                else:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep + (timestep ** 2) * force / 2)

            # force after position update (needed for verlet velocity)
            new_rel_pos, new_rel_dis = atomic_distances(pos_steps[i + 1, :, :], box_dim)
            new_force = lj_force(new_rel_pos, new_rel_dis)[1]

            if dimless:
                vel_steps[i + 1, :, :] = vel_steps[i, :, :] + timestep * (new_force + force) / 2
            else:
                vel_steps[i + 1, :, :] = vel_steps[i, :, :] + timestep * (new_force + force) / 2 / ARG_MASS
        else:
            # Euler
            if periodic:
                # make sure it's inside box dimension -> modulus gives periodicity
                if dimless:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep) % box_dim
                else:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep) % box_dim
            else:
                if dimless:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep)
                else:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep)

            rel_pos = atomic_distances(pos, box_dim)[0]
            rel_dis = atomic_distances(pos, box_dim)[1]
            force = lj_force(rel_pos, rel_dis)[1]

            if dimless:
                vel_steps[i + 1, :, :] = vel_steps[i, :, :] + force * timestep
            else:
                vel_steps[i + 1, :, :] = vel_steps[i, :, :] + force * timestep / ARG_MASS

        if rescaling and (i %rescaling_timesteps==0) and (i < rescaling_max_timesteps+1):
            # Rescale velocity
            if rescaling_mode == 0:
                # old kin energy avg
                rescaling1 = np.sum([kinetic_energy(vel_steps[i-x, :, :])[1] for x in range(min(i+1,5000))])/min(i+1,5000)
                # new kin energy avg
                rescaling2 = np.sum([kinetic_energy(vel_steps[i+1-x, :, :])[1] for x in range(min(i + 1, 5000))]) / min(i + 1, 5000)
                # rescaling factor (in sqrt(...) so values get closer to 1)
                v_lambda = np.sqrt((num_atoms - 1) * 3 * KB * temp / (EPSILON * np.sum([np.sqrt(np.sum([v[i] ** 2 for i in range(dim)])) for v in vel_steps[i + 1, :, :]])))  # / TEMP
                #v_lambda = np.sqrt((num_atoms - 1) * 3 * KB * temp / (ARG_MASS * np.sum([np.sqrt(np.sum([v[i] ** 2 for i in range(dim)])) for v in vel_steps[i + 1, :, :]]) * dimless_velocity))/1500
                #v_lambda = np.sqrt((num_atoms - 1) * 3 * KB * temp / (ARG_MASS * np.sum([np.sqrt(np.sum([v[i] ** 2 for i in range(dim)])) for v in vel_steps[i + 1, :, :]]))) / TEMP
                current_temperature = rescaling1 * EPSILON / ((num_atoms - 1) * 3 / 2 * KB)
                need_rescaling = np.abs(rescaling2 - rescaling1) < rescaling_delta*0.015
            else:
                # target kin energy
                rescaling1 = (num_atoms - 1) * 3 / 2 * temp * KB / EPSILON
                # new kin energy avg
                rescaling2 = np.sum([kinetic_energy(vel_steps[i+1-x, :, :])[1] for x in range(min(i + 1, 5000))]) / min(i + 1, 5000)
                # rescaling factor (in sqrt(...) so values get closer to 1)
                v_lambda = np.sqrt(np.sqrt((num_atoms - 1) * 3 / 2 * KB * temp / (EPSILON * np.sum([np.sqrt(np.sum([v[i] ** 2 for i in range(dim)])) for v in vel_steps[i + 1, :, :]]))))  # / TEMP
                current_temperature = rescaling2 * EPSILON / ((num_atoms - 1) * 3 / 2 * KB)
                need_rescaling = np.abs(rescaling2 - rescaling1) > rescaling_delta * current_temperature

            if need_rescaling:
                # limit rescaling factor between 0.5 and 2.0
                v_lambda = max(0.5,v_lambda)
                v_lambda = min(2.0,v_lambda)
                # apply rescaling factor
                vel_steps[i + 1, :, :] *= v_lambda
                # rescaling statistics below
                rescale_counter+=1
                rescale_max = max(rescale_max,v_lambda)
                rescale_min = min(rescale_min,v_lambda)

    if rescaling:
        # print rescaling statistics
        print("Rescaled",rescale_counter,"times with λ: [",rescale_min,"~",rescale_max,"]")

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
                        if dis > (box_dim * 0.5):
                            dis = dis - box_dim
                        if dis <= -(box_dim * 0.5):
                            dis = dis + box_dim
                    else:
                        if dis > (box_dim * 0.5):
                            dis = dis - box_dim
                        if dis <= -(box_dim * 0.5):
                            dis = dis + box_dim
                rel_pos[i][j][k] = dis

    rel_dist = np.zeros([len(pos), len(pos)])
    for i in range(0, len(rel_pos)):
        for j in range(0, len(rel_pos)):
            rel_dist[i][j] = np.math.sqrt(sum(i ** 2 for i in rel_pos[i][j]))

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

    if dimless:
        for i in range(0, len(rel_pos[1])):  # particle i
            for j in range(0, len(rel_pos[1])):  # particle i rel to j (!=i)
                if i != j:
                    dUdt[i, j] = -24 * ((2 / (rel_dist[i, j] ** 13)) - (1 / rel_dist[i, j] ** 7)) / (
                                     rel_dist[i, j])
                else:
                    dUdt[i, j] = 0
        for i in range(0, len(rel_pos[1])):  # particle i
            for j in range(0, len(rel_pos[1])):  # particle i rel to j (!=i)
                force[i, j, :] = dUdt[i, j] * rel_pos[i, j, :]
                
    else:
        for i in range(0, len(rel_pos[1])):  # particle i
            for j in range(0, len(rel_pos[1])):  # particle i rel to j (!=i)
                if i != j:
                    dUdt[i, j] = -24 * EPSILON * (
                            (2 * SIGMA ** 12 / (rel_dist[i, j] ** 13)) - (SIGMA ** 6 / rel_dist[i, j] ** 7)) / (
                                     rel_dist[i, j])
                else:
                    dUdt[i, j] = 0
        for i in range(0, len(rel_pos[1])):  # particle i
            for j in range(0, len(rel_pos[1])):  # particle i rel to j (!=i)
                force[i, j, :] = dUdt[i, j] * rel_pos[i, j, :]

    # while this looks horrible, and is horrible, it works. However, needs significant optimazation
        
    force_atom = np.sum(force, axis=1)

    return force, force_atom


def lj_force_old(rel_pos, rel_dist):
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

    for i in range(0, len(rel_pos[1])):  # particle i
        for j in range(0, len(rel_pos[1])):  # particle i rel to j (!=i)
            if i != j:
                dUdt[i, j] = -24 * EPSILON * (
                        (2 * SIGMA ** 12 / (rel_dist[i, j] ** 13)) - (SIGMA ** 6 / rel_dist[i, j] ** 13)) / (
                                 rel_dist[i, j])
            else:
                dUdt[i, j] = 0
    for i in range(0, len(rel_pos[1])):  # particle i
        for j in range(0, len(rel_pos[1])):  # particle i rel to j (!=i)
            force[i, j, :] = dUdt[i, j] * rel_pos[i, j, :]
    # while this looks horrible, and is horrible, it works. However, needs significant optimazation

    force_atom = np.sum(force, axis=1)

    return force, force_atom


def fcc_lattice_old(num_atoms, lat_const):
    """
    Initializes a system of atoms on an fcc lattice.
    
    NOTE CURRENTLY, ONLY WORKS FOR 4 ATOMS
    Initial vectors are:
    a1 = [D,0,0]
    a2 = [0,D,0]
    a3 = [0,0,D]
    Here, D is the distance between 2 adjecent corner atoms.
    
    lattice basis vectors are:
    r1 = [0,0,0]
    r2 = 1/2(a1+a2)
    r3 = 1/2(a2+a3)
    r4 = 1/2(a3+a1)
    
    FCC lattice is only possible in 3D due to definition of FCC lattice
    
    https://solidstate.quantumtinkerer.tudelft.nl/10_xray/ can be used as a reference

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
    if num_atoms % 4 == 0:
        a = np.array([[lat_const, 0, 0], [0, lat_const, 0], [0, 0, lat_const]])
        BZ = np.int(num_atoms/4)
        print('N = multiple of 4')
        # below is not elegant at all, but it works without writing over complex code for a simple thing.
        pos_vec = 0.5 * np.array([[0., 0., 0.], np.add(a[0, :], a[1, :]), np.add(a[1, :], a[2, :]), np.add(a[2, :], a[0, :])])
        # offset can be usefull for plotting purposes. Update required to match boxsize regarding offset
        offset = [0, 0, 0] 
        pos_vec = np.add(pos_vec, offset)
        print('fcc lattice vector is', pos_vec)
    else:
        print('N is not multiple of 4, FCC lattice not possible ')
        exit()
    return pos_vec


def fcc_lattice(num_atoms, lat_const):
    """
    Initializes a system of atoms on an fcc lattice.
    
    NOTE CURRENTLY, ONLY WORKS FOR 4 ATOMS
    Initial vectors are:
    a1 = [D,0,0]
    a2 = [0,D,0]
    a3 = [0,0,D]
    Here, D is the distance between 2 adjecent corner atoms.
    
    lattice basis vectors are:
    r1 = [0,0,0]
    r2 = 1/2(a1+a2)
    r3 = 1/2(a2+a3)
    r4 = 1/2(a3+a1)
    
    FCC lattice is only possible in 3D due to definition of FCC lattice
    
    https://solidstate.quantumtinkerer.tudelft.nl/10_xray/ can be used as a reference

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
    if num_atoms >=4:
        a = np.array([[lat_const, 0, 0], [0, lat_const, 0], [0, 0, lat_const]])
        BZ = np.int(num_atoms/4)
        print('N = multiple of 4')
        # below is not elegant at all, but it works without writing over complex code for a simple thing.
        pos_vec = 0.5 * np.array([[0., 0., 0.], np.add(a[0, :], a[1, :]), np.add(a[1, :], a[2, :]), np.add(a[2, :], a[0, :])])
        # offset can be usefull for plotting purposes. Update required to match boxsize regarding offset
        offset = [0, 0, 0] 
        pos_vec = np.add(pos_vec, offset)
        print(pos_vec)
        print(a[0,:])
        if num_atoms>4:
            for i in range(2):
                pos_ext = pos_vec+a[i,:]
                pos_vec = np.append(pos_vec, pos_ext, axis=0)
            pos_vec = np.append(pos_vec, pos_vec+a[2,:], axis=0)
            print('fcc lattice vector is', pos_vec)
            
    else:
        print('N is not multiple of 4, FCC lattice not possible ')
        exit()
    return pos_vec


vector = fcc_lattice(32, 1)
print(vector.shape)


def locationplot(locations, latmult):
    """
    Plots locations of N particles
    

    Parameters
    ----------
    locations : np.ndarray
        locations of particles
    latmult: scalar
        How many lattices are to be plotted\
        latmult=1 4particles; latmult=2 16particles
    
    Returns
    -------
    plot : plt.plot
        plot of the locations of the particles.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(locations[:, 0], locations[:, 1], locations[:, 2])
    ax.set_xlim3d(0, latmult)
    ax.set_ylim3d(0, latmult)
    ax.set_zlim3d(0, latmult)
    plt.show()



locationplot(vector,2)


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
            velsquared = vel**2 ##np.power(vel, 2.0)
            vel_summed = np.sum(velsquared, axis=1)
            vel_abs = vel_summed**0.5 #np.power(vel_summed, 0.5) #the total velocity, of 1 particle stored in an array for each particle. Since a bug was present, This is rewritten in, over simplified steps.
            ke_part = 0.5 * vel_abs**2 #np.power(vel_abs, 2)
            ke_total = np.sum(ke_part)

        else:
            ke += 0.5 * ARG_MASS * np.power(np.math.sqrt(sum(i ** 2 for i in vel[i])), 2.0)

    return ke_part, ke_total


def kinetic_energyold(vel):
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
            ke += 0.5 * np.power(np.math.sqrt(sum(i ** 2 for i in vel[i])), 2.0)
        else:
            ke += 0.5 * ARG_MASS * np.power(np.math.sqrt(sum(i ** 2 for i in vel[i])), 2.0)

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
                pot_e[i][j] = 4 * EPSILON * ((SIGMA / rel_dist[i][j]) ** 12 - (SIGMA / rel_dist[i][j]) ** 6)
            else:
                pot_e[i][j] = 0
    pot_e_particle = np.sum(pot_e, axis=1)
    pot_total = np.sum(pot_e_particle) / 2

    if dimless:
        for j in range(0, num_atoms1):
            for i in range(0, num_atoms1):
                if i != j:
                    pot_e[i][j] = 4 * ((1 / rel_dist[i][j]) ** 12 - (1 / rel_dist[i][j]) ** 6)
                else:
                    pot_e[i][j] = 0
        pot_e_particle = np.sum(pot_e, axis=1)
        pot_total = np.sum(pot_e_particle) / 2

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

    return kinetic_energy(vel)[1] + potential_energy(rel_dist)[2]


def ms_displacement(loc, timestep):
    """
        Computes the mean square displacement of a single atom.

        Parameters
        ----------
        loc: np.ndarray
            locations of particles over time [timestep, particle, dims]
        timestep : int
            the timestep of the particle which is used as initial value

        Returns
        -------
        msd_1: np.ndarray
            The msq time dependent array, for N dimensions and M particles 
        msd_2: np.ndarray
            the msq time dependent array, summed over the dimensions, for M particles
            [msd_part1(dtime=0), msd_part2(dtime=0),... ], [msd_part1(dtime=1), msd_part2(time=1),... ], ....
        msd_3: np.ndarray
            the msq time dependent vector, summed over both dimensions and particles
            [msd_total(dtime=0), msd_total(dtime=1),.....]
        -------

        """
    init_loc = loc[timestep, :, :]
    loc_usage = loc[timestep:-1, :, :]
    msd_1 = np.abs((loc_usage - init_loc)**2)
    # next
    for i in range(len(loc_usage[:,0,0])):
        msd_1[i, :, :] = msd_1[i,:,:] / (i + 1)
        
    msd_2 = np.sum(msd_1, axis=2)
    N = len(program[0][0,:,0]) #number of particles
    msd_3 = np.sum(msd_2, axis=1)/N
    print(len(msd_3))
    D = np.empty(len(msd_3))
    for i in range(len(msd_3)):
        D[i] = msd_3[i] / (6*(i + 1))
        
    plt.plot(D)
    plt.show()
    print('the diff coeff is shown in the plot', D)
    return msd_1, msd_2, msd_3


def msd_plot(msd,partnum):
    """"
    plots the MSD of a single atom NOTE MIGHT NOW WORK
    Parameters
    ----------
    msd_2: np.ndarray
        the msq time dependent array, summed over the dimensions, for M particles
        [msd_part1(dtime=0), msd_part2(dtime=0),... ], [msd_part1(dtime=1), msd_part2(time=1),... ], ....
        best use case: msd_2 from ms_displacement function
    partnum: int
        the particle that is to be plotted by the function
    Returns
    ----------
    None
    """
    plt.plot(msd[:,0])
    plt.show()
    return


def process_data():
    print("Test if the total energy is conserved")
    pos1 = positions_store[0, :, :]
    pos2 = positions_store[steps - 1, :, :]

    vel1 = velocities_store[0, :, :]
    vel2 = velocities_store[steps - 1, :, :]

    r_pos1 = atomic_distances(pos1, box_dim)
    r_pos2 = atomic_distances(pos2, box_dim)

    print("Initial total energy: " + str(total_energy(vel1, r_pos1[1])))
    print("Final total energy:   " + str(total_energy(vel2, r_pos2[1])))
    print("Delta total energy:   " + str(total_energy(vel2, r_pos2[1]) - total_energy(vel1, r_pos1[1])))

    times = np.linspace(0, dt * steps, steps)
    if num_atoms == 2:
        print("Plot inter-atom distance over time")
        if dimless:
            distances = [np.max(atomic_distances(positions_store[x, :, :], box_dim)[1]) for x in range(steps)]
        else:
            distances = [np.max(atomic_distances(positions_store[x, :, :], box_dim)[1]) for x in range(steps)]
        plt.plot(times, distances)
        if dimless:
            plt.ylabel('Distance (dimless)')
            plt.xlabel('Time (dimless)')

        else:
            plt.ylabel('Distance (m)')
            plt.xlabel('Time (s)')

        plt.show()

    print("Print energy levels over time")
    if dimless:
        energies = [(kinetic_energy(velocities_store[x, :, :])[1],
                     potential_energy(atomic_distances(positions_store[x, :, :], box_dim)[1])[2],
                     total_energy(velocities_store[x, :, :],
                                  atomic_distances(positions_store[x, :, :], box_dim)[1]))
                    for x in range(steps)]
    else:
        energies = [kinetic_energy(velocities_store[x, :, :])[1] for x in range(steps)]
    # times = np.linspace(0, dt*steps, steps)
    plt.plot(times, energies)
    plt.xlabel('Time (dimless)')
    plt.ylabel('Energy (dimless)')
    plt.legend(('kinetic energy', 'potential energy', 'total energy'))
    plt.show()

def test_initial_velocities(init_velocities):
    init_velocities = init_velocity(1000, TEMP, dim)

    vel_mag = [np.sqrt(np.sum([v[i] ** 2 for i in range(dim)])) for v in init_velocities]

    plt.hist(vel_mag, bins=15)
    plt.show()

def main():
    """"
        Beginning of program
        to start with "easy" locations and velocities,
        init_pos = [[0, 0], [0, 1]]
        init_vel = [[1, 1], [1, 1]]

    """
    #    easy, handpicked initial positions and velocities.
    # init_pos = [[9.9, 9.6, 8.7], [0.3, 0.6, 0.3], [3.5, 4.6, 5.7], [9.9, 3.3, 6.6], [6.0, 7.5, 9.0],
    #            [0.6, 0.6, 9.0], [3.3, 3.3, 3.3], [8.8, 2.7, 6.3], [6.3, 8.7, 1.5]]
    # init_vel = [[1.2, 0.8, 1.2], [-0.9, -0.67, -0.88], [-0.89, 0.94, 1.55], [1.52, -0.53, 0.97], [0.60, -1.55, 0.32]
    #    , [-0.22, 1.53, -0.34], [1.25, 0.66, -0.97], [-0.36, -1.29, 0.09], [1.22, 0.01, -0.61]]

    # Below is the must be uncommented for the delivarble.
    init_pos = fcc_lattice(num_atoms, 1.547)
    init_vel = init_velocity(num_atoms,TEMP,dim)

    #init_pos = [[25,25,25], [28,25,25], [25,25,27]]
    #init_vel = [[1,0,0], [1,0,0], [1,0,0]]

    test_initial_velocities(init_vel)

    simulate(init_pos, init_vel, steps, dt, box_dim)
    process_data()
    p1 = positions_store
    v1 = velocities_store


    global verlet
    verlet = False

    simulate(init_pos, init_vel, steps, dt, box_dim)
    process_data()
    p2 = positions_store
    v2 = velocities_store

    print("Maximum error in positions data: ", np.max(p2-p1))
    print("Maximum error in velocities data: ", np.max(v2-v1))
    return p1, v1, p2, v2

program = main()

# original positions
plt.plot(program[0][:,0,0])
plt.show()

# make positions continuous
p0 = (box_dim/2 - np.abs(box_dim/2-program[0][:,0,0]))
plt.plot(p0)
plt.show()

qq = ms_displacement(program[0], 100)
plt.plot(qq[0][:,0,0]) #plot of msd particle 0 in 0 axis
plt.show()

plt.plot(program[0][:,0,0]) #location of particle 0 in 0 axis
plt.show()

a = ms_displacement(program[0], 10000)

plt.plot(a[2])


