# coding: utf-8

"""
This is a suggestion for structuring your simulation code properly.
However, it is not set in stone. You may modify it if you feel like
you have a good reason to do so.
"""

import numpy as np
import matplotlib.pyplot as plt
import time

# initalizing self defined system parameters
num_atoms = 4  # amount of particles
dim = 3  # dimensions
box_dim = 3.313  # 2* 1.547 #10  # meters; bounding box dimension
dt = 1e-4  # s; stepsize

steps = 10000  # amount of steps
dimless = True  # use dimensionless units
periodic = True  # use periodicity
verlet = True  # use Verlet's algorithm (false: Euler's algorithm)
rescaling = True  # use Temperature rescaling
rescaling_mode = 1  # 0 = kin-NRG-based | temp-based
rescaling_delta = 0.0027  # delta for activation of rescaling
rescaling_timesteps = steps / 10  # timesteps interval for rescaling check
rescaling_max_timesteps = steps / 2  # max timesteps for rescaling
rescaling_limit = True  # rescale limit [lower~upper]
rescaling_limit_lower = 0.5
rescaling_limit_upper = 2.0
rescaling_factor = 0.5  # 0.5 = sqrt

# +
# for 32 particles
# Lattice constant = box_dim/2
# liquid
# box_dim = 3.313
# T = 1

# solid
# box_dim = 3.41995
# T = 0.5

# Gas
# 4.7425
# T = 3

# -

# Parameters physical, supplied by course, or related to Argon
temp = 2  # Kelvin
TEMP = 119.8  # Kelvin
KB = 1.38064852e-23  # m^2*kg/s^2/K
SIGMA = 3.405e-10  # meter
EPSILON = TEMP * KB  # depth of potential well/dispersion energy
N_b = 6.02214076e23  # Avogadros number; 1/mol
R = 8.31446261815324  # J/K/mole; universal gas constant
ARG_UMASS = 39.95  # u; atomic mass of argon
ARG_MMASS = ARG_UMASS / 1000  # kg/mol; mole mass of argon
ARG_MASS = ARG_UMASS * 1.6605e-27  # Kg mass of a single atom in Kg

# conversion values for dimensionless units
dimless_time = 1.0 / np.math.sqrt((ARG_MASS * SIGMA ** 2 / EPSILON))  # s; dimensionless time
dimless_energy = 1.0 / EPSILON  # J; dimensionless energy
dimless_distance = 1.0 / SIGMA  # m; dimensionless distance
dimless_velocity = 1.0 / np.math.sqrt(EPSILON / ARG_MASS)  # m/s; dimensionless velocity


def init_velocity(number_atoms, temperature, dimensions):
    """
    Initializes the system with Gaussian distributed velocities. This
    init_velocity is loosely based on 3D system, however it will output
    2D just fine, although more pertubated. This function is based a
    simplified Boltzmann distribution, found at:
    https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann
    _distribution#Typical_speeds
    Parameters
    ----------
    number_atoms : int
        The number of particles in the system.
    temperature : float
        The (unitless) temperature of the system.
    dimensions : int
        The dimensions of the system.

    Returns
    -------
    vel_vec : np.ndarray
        Array of particle velocities
    """
    vel_p = np.sqrt(2 * KB * temperature / ARG_MASS)

    if dimless:
        vel_p *= dimless_velocity

    vel_mean = 2 * vel_p / np.sqrt(np.pi)
    vel_msq = (3 * vel_p ** 2) / 2
    vel_std = vel_msq - (vel_mean ** 2)
    vel_vec = np.random.normal(vel_mean, vel_std, (number_atoms, dimensions))
    vel_mag = np.linalg.norm(vel_vec, axis=1)
    vel_vec *= vel_mean / np.mean(vel_mag)  # Rescale the magnitudes to match the vel_mean speed

    # create random negativity
    for v in range(number_atoms):
        for i in range(dimensions):
            # either *1 or *-1
            vel_vec[v, i] *= (1 - 2 * np.random.randint(2))

    vel_vec -= np.mean(vel_vec)  # remove mean for no drift velocity

    return vel_vec


def init_position(number_atoms, box_dimensions, dimensions):
    """
    Initializes the system with random positions.
    This does not require non dimensionalization scaling, since it is
    not based on physical parameters.

    Parameters
    ----------
    number_atoms : int
        The number of particles in the system.
    box_dimensions : float
        The dimension of the simulation box
    dimensions : int
        The dimensions of the system.

    Returns
    -------
    pos_vec : np.ndarray
        Array of particle positions
    """
    randoms = np.random.random((number_atoms, dimensions))
    pos_vec = randoms * box_dimensions

    return pos_vec


def simulate(init_pos, init_vel, num_tsteps, timestep, box_dimensions):
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
    box_dimensions : float
        Dimensions of the simulation box

    Returns
    -------
    Any quantities or observables that you wish to study.
    """

    if verlet:
        print("Simulating using Verlet's algorithm")
    else:
        print("Simulating using Euler's algorithm")

    pos_steps = np.zeros((num_tsteps, num_atoms, dim))
    vel_steps = np.zeros((num_tsteps, num_atoms, dim))

    pos_steps[0, :, :] = init_pos
    vel_steps[0, :, :] = init_vel

    rescale_counter = 0
    rescale_max = 1.0  # WHY DIT?
    rescale_min = 1.0

    for i in range(num_tsteps - 1):
        pos = pos_steps[i, :, :]

        if verlet:
            rel_pos, rel_dis = atomic_distances(pos, box_dimensions)
            force = lj_force(rel_pos, rel_dis)[1]

            # Keep particle inside box using modulus when periodic
            if periodic:
                if dimless:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep
                                              + (timestep ** 2) * force / 2) % box_dimensions
                else:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep + (
                            timestep ** 2) * force / 2) % box_dimensions

            else:
                if dimless:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep + (timestep ** 2) * force / 2)
                else:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep + (timestep ** 2) * force / 2)

            # force after position update (needed for verlet velocity)
            new_rel_pos, new_rel_dis = atomic_distances(pos_steps[i + 1, :, :], box_dimensions)
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
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep) % box_dimensions
                else:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep) % box_dimensions
            else:
                if dimless:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep)
                else:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep)

            rel_pos = atomic_distances(pos, box_dimensions)[0]
            rel_dis = atomic_distances(pos, box_dimensions)[1]
            force = lj_force(rel_pos, rel_dis)[1]

            if dimless:
                vel_steps[i + 1, :, :] = vel_steps[i, :, :] + force * timestep
            else:
                vel_steps[i + 1, :, :] = vel_steps[i, :, :] + force * timestep / ARG_MASS

        if rescaling and (i % rescaling_timesteps == 0) and (i < rescaling_max_timesteps + 1):
            # Rescale velocity
            if rescaling_mode == 0:
                # old kin energy avg
                rescaling1 = np.sum([kinetic_energy(vel_steps[i - x, :, :])[1] for x in range(min(i + 1, 5000))]) / min(
                    i + 1, 5000)
                # new kin energy avg
                rescaling2 = np.sum(
                    [kinetic_energy(vel_steps[i + 1 - x, :, :])[1] for x in range(min(i + 1, 5000))]) / min(i + 1, 5000)
                # rescaling factor (in sqrt(...) so values get closer to 1)
                v_lambda = np.sqrt((num_atoms - 1) * 3 * KB * temp / (EPSILON * np.sum(
                    [np.sqrt(np.sum([v[i] ** 2 for i in range(dim)])) for v in vel_steps[i + 1, :, :]])))  # / TEMP
                # v_lambda = np.sqrt((num_atoms - 1) * 3 * KB * temp / (ARG_MASS * np.sum([np.sqrt(np.sum([v[i] ** 2
                # for i in range(dim)])) for v in vel_steps[i + 1, :, :]]) * dimless_velocity))/1500
                # v_lambda = np.sqrt((num_atoms - 1) * 3 * KB * temp / (ARG_MASS * np.sum([np.sqrt(np.sum([v[i] ** 2
                # for i in range(dim)])) for v in vel_steps[i + 1, :, :]]))) / TEMP
                # current_temperature = rescaling1 * EPSILON / ((num_atoms - 1) * 3 / 2 * KB)
                need_rescaling = np.abs(rescaling2 - rescaling1) < rescaling_delta * 0.015
            else:
                # target kin energy
                rescaling1 = (num_atoms - 1) * 3 / 2 * temp * KB / EPSILON
                # new kin energy avg
                kin_nrg = np.zeros(int(min(i + 1, int(rescaling_timesteps))))
                for x in range(len(kin_nrg)):
                    kin_nrg[x] = kinetic_energy(vel_steps[i + 1 - x, :, :])[1]
                rescaling2 = np.sum(kin_nrg) / int(min(i + 1, int(rescaling_timesteps)))
                # rescaling factor (in sqrt(...) so values get closer to 1)
                # v_lambda = np.sqrt(np.sqrt((num_atoms - 1) * 3 / 2 * KB * temp / (EPSILON * np.sum(
                # [np.sqrt(np.sum([v[i] ** 2 for i in range(dim)])) for v in vel_steps[i + 1, :, :]]))))  # / TEMP
                v_lambda = np.power((num_atoms - 1) * 3 / 2 * KB * temp / (
                        EPSILON * np.sum(np.linalg.norm(vel_steps[i + 1, :, :], axis=1))),
                                    rescaling_factor)  # / TEMP
                current_temperature = rescaling2 * EPSILON / ((num_atoms - 1) * 3 / 2 * KB)
                need_rescaling = np.abs(rescaling2 - rescaling1) > rescaling_delta * current_temperature

            if need_rescaling:
                # limit rescaling factor between 0.5 and 2.0
                if rescaling_limit:
                    v_lambda = max(rescaling_limit_lower, v_lambda)
                    v_lambda = min(rescaling_limit_upper, v_lambda)

                # apply rescaling factor
                vel_steps[i + 1, :, :] *= v_lambda
                # rescaling statistics below
                rescale_counter += 1
                rescale_max = max(rescale_max, v_lambda)
                rescale_min = min(rescale_min, v_lambda)

    if rescaling:
        # print rescaling statistics
        print("Rescaled", rescale_counter, "times with λ: [", rescale_min, "~", rescale_max, "]")

    return pos_steps, vel_steps


def atomic_distances(pos, box_dimensions):
    """
    Calculates relative positions and distances between particles.

    parameters
    ----------
    pos : np.ndarray
        The positions of the particles in cartesian space
    box_dimensions : float
        The dimension of the simulation box

    returns
    -------
    rel_pos : np.ndarray
        Relative positions of particles
    rel_dist : np.ndarray
        The distance between particles
    """

    dimensions = len(pos[0])
    # NOTE: includes rel_dist/rel_pos to itself (= 0 / [0.0, 0.0])
    rel_pos = np.zeros([len(pos), len(pos), dimensions])

    for i in range(0, len(pos)):
        for j in range(0, len(pos)):
            for k in range(0, dimensions):
                dis = pos[j][k] - pos[i][k]
                if periodic:
                    if dimless:
                        if dis > (box_dimensions * 0.5):
                            dis = dis - box_dimensions
                        if dis <= -(box_dimensions * 0.5):
                            dis = dis + box_dimensions
                    else:
                        if dis > (box_dimensions * 0.5):
                            dis = dis - box_dimensions
                        if dis <= -(box_dimensions * 0.5):
                            dis = dis + box_dimensions
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
    dudt = np.zeros([len(rel_dist), len(rel_dist)])
    force = np.zeros([len(rel_pos[1]), len(rel_pos[1]), len(rel_pos[0][0])])

    if dimless:
        for i in range(0, len(rel_pos[1])):  # particle i
            for j in range(0, len(rel_pos[1])):  # particle i rel to j (!=i)
                if i != j:
                    dudt[i, j] = -24 * ((2 / (rel_dist[i, j] ** 13)) - (1 / rel_dist[i, j] ** 7)) / (
                        rel_dist[i, j])
                else:
                    dudt[i, j] = 0
        for i in range(0, len(rel_pos[1])):  # particle i
            for j in range(0, len(rel_pos[1])):  # particle i rel to j (!=i)
                force[i, j, :] = dudt[i, j] * rel_pos[i, j, :]

    else:
        for i in range(0, len(rel_pos[1])):  # particle i
            for j in range(0, len(rel_pos[1])):  # particle i rel to j (!=i)
                if i != j:
                    dudt[i, j] = -24 * EPSILON * (
                            (2 * SIGMA ** 12 / (rel_dist[i, j] ** 13)) - (SIGMA ** 6 / rel_dist[i, j] ** 7)) / (
                                     rel_dist[i, j])
                else:
                    dudt[i, j] = 0
        for i in range(0, len(rel_pos[1])):  # particle i
            for j in range(0, len(rel_pos[1])):  # particle i rel to j (!=i)
                force[i, j, :] = dudt[i, j] * rel_pos[i, j, :]

    # while this looks horrible, and is horrible, it works. However, needs significant optimazation

    force_atom = np.sum(force, axis=1)

    return force, force_atom


def fcc_lattice(number_atoms, lat_const):
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
    number_atoms : int
        The number of particles in the system
    lat_const : float
        The lattice constant for an fcc lattice

    Returns
    -------
    pos_vec : np.ndarray
        Array of particle coordinates
    """
    # placeholder
    pos_vec = 0

    if number_atoms >= 4:
        a = np.array([[lat_const, 0, 0], [0, lat_const, 0], [0, 0, lat_const]])
        # BZ = int(number_atoms / 4)
        print('N = multiple of 4')
        # below is not elegant at all, but it works without writing over complex code for a simple thing.
        pos_vec = 0.5 * np.array(
            [[0., 0., 0.], np.add(a[0, :], a[1, :]), np.add(a[1, :], a[2, :]), np.add(a[2, :], a[0, :])])
        # offset can be usefull for plotting purposes. Update required to match boxsize regarding offset
        offset = [0.01 * box_dim, 0.01 * box_dim, 0.01 * box_dim]  # NOTE I ADDED OFFSET
        pos_vec = np.add(pos_vec, offset)
        print(pos_vec)
        # print(a[0,:])
        if number_atoms > 4:
            for i in range(2):
                pos_ext = pos_vec + a[i, :]
                pos_vec = np.append(pos_vec, pos_ext, axis=0)
            pos_vec = np.append(pos_vec, pos_vec + a[2, :], axis=0)
            print('fcc lattice vector is', pos_vec)

    else:
        print('N is not multiple of 4, FCC lattice not possible ')
        exit()
    return pos_vec


def fcc_lattice_big(number_atoms, lat_const):
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
    number_atoms : int
        The number of particles in the system
    lat_const : float
        The lattice constant for an fcc lattice

    Returns
    -------
    pos_vec : np.ndarray
        Array of particle coordinates
    """
    # placeholder
    pos_vec = 0
    if number_atoms == 14:
        a = np.array([[lat_const, 0, 0], [0, lat_const, 0], [0, 0, lat_const]])
        # BZ = int(number_atoms / 4)
        print('N = multiple of 4')
        # below is not elegant at all, but it works without writing over complex code for a simple thing.
        pos_vec = 0.5 * np.array(
            [[0., 0., 0.], np.add(a[0, :], a[1, :]), np.add(a[1, :], a[2, :]), np.add(a[2, :], a[0, :])])
        pos_vec_ext = np.array([a[0, :], a[1, :], a[2, :], a[0, :] + a[1, :], a[0, :] + a[2, :], a[1, :] + a[2, :],
                                a[0, :] + a[1, :] + a[2, :], a[0, :] + 0.5 * (a[1, :] + a[2, :]),
                                a[1, :] + 0.5 * (a[2, :] + a[0, :]), a[2, :] + 0.5 * (a[1, :] + a[0, :])])
        pos_vec = np.append(pos_vec, pos_vec_ext, axis=0)
    else:
        print('value not 14')
    print('fcc lattice vector is', pos_vec)
    return pos_vec


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


# q = fcc_lattice(32,1)
# locationplot(q,2)


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

    if dimless:
        velsquared = vel ** 2  # np.power(vel, 2.0)
        vel_summed = np.sum(velsquared, axis=1)
        vel_abs = vel_summed ** 0.5  # np.power(vel_summed, 0.5)
        # the total velocity, of 1 particle stored in an array for each particle.
        # Since a bug was present, This is rewritten in, over simplified steps.
        ke_part = 0.5 * vel_abs ** 2  # np.power(vel_abs, 2)
        ke_total = np.sum(ke_part)
    else:
        ke = 0
        for i in range(0, len(vel)):
            ke += 0.5 * ARG_MASS * np.power(np.math.sqrt(sum(i ** 2 for i in vel[i])), 2.0)
        return ke, ke

    return ke_part, ke_total


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

    return (kinetic_energy(vel))[1] + (potential_energy(rel_dist))[2]