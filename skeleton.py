# coding: utf-8
"""
This is a suggestion for structuring your simulation code properly.
However, it is not set in stone. You may modify it if you feel like
you have a good reason to do so.
"""

# +
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import time
# LOAD SEED

from pickle import load
# -

start_time = time.time()

# initalizing self defined system parameters
num_atoms = 32  # amount of particles
dim = 3  # dimensions
box_dim = 3.313  # 2* 1.547 #10  # meters; bounding box dimension
dt = 1e-4  # s; stepsize

steps = 10000  # amount of steps
dimless = True  # use dimensionless units
periodic = True  # use periodicity
verlet = True  # use Verlet's algorithm (false: Euler's algorithm)
rescaling = True  # use Temperature rescaling
rescaling_mode = 1  # 0 = kin-NRG-based | temp-based
rescaling_delta_energy = 0.0000001 * num_atoms # delta for activation of rescaling
rescaling_max_timesteps = 5000  # max timesteps for rescaling
rescaling_max_rescales = 60
rescaling_limit = True  # rescale limit [lower~upper]
rescaling_limit_lower = 0.5
rescaling_limit_upper = 2.0
rescaling_factor = 0.5  # 0.5 = sqrt
average_rescale = 100
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
temp = 1  # Kelvin
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

with open('state.obj', 'rb') as f:
    np.random.set_state(load(f))


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

    pos_steps = np.zeros((num_tsteps, num_atoms, dim))
    vel_steps = np.zeros((num_tsteps, num_atoms, dim))
    kin_steps = np.zeros(num_tsteps)
    energy_steps = np.zeros((3, num_tsteps))
    
    pos_steps[0, :, :] = init_pos
    vel_steps[0, :, :] = init_vel
    kin_steps[0] = kinetic_energy(init_vel)[1]
    
    rel_pos, rel_dis = atomic_distances(init_pos, box_dim)
    energy_steps[0, 0] = kinetic_energy(init_vel)[1]
    energy_steps[1, 0] = potential_energy(rel_dis)[2]
    energy_steps[2, 0] = energy_steps[0, 0]+energy_steps[1, 0]
    
    init_rel_dist = atomic_distances(init_pos, box_dim)
    print("Initial total energy: " + str(total_energy(init_vel, init_rel_dist[1])))
    
    
 

    
    rescale_counter = 0
    rescale_max = 1.0                                          #WHY DIT?
    rescale_min = 1.0
    rescaling_moment = 10


    for i in range(num_tsteps-1):
        pos = pos_steps[i, :, :]

        if verlet:
            rel_pos, rel_dis = atomic_distances(pos, box_dim)
            force = lj_force(rel_pos, rel_dis)[1]

            # Keep particle inside box using modulus when periodic
            if periodic:
                if dimless:
                    pos_steps[i+1, :, :] = (pos + vel_steps[i, :, :]*timestep
                                            + (timestep**2) * force / 2) % box_dim
                else:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep + (
                            timestep ** 2) * force / 2) % box_dim

            else:
                if dimless:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep + (timestep ** 2) * force / 2)
                else:
                    pos_steps[i + 1, :, :] = (pos + vel_steps[i, :, :] * timestep + (timestep ** 2) * force / 2)

            # force after position update (needed for verlet velocity)
            rel_pos, rel_dis = atomic_distances(pos_steps[i + 1, :, :], box_dim)
            new_force = lj_force(rel_pos, rel_dis)[1]

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

            rel_pos = atomic_distances(pos_steps[i + 1, :, :], box_dim)[0]
            rel_dis = atomic_distances(pos_steps[i + 1, :, :], box_dim)[1]
            force = lj_force(rel_pos, rel_dis)[1]

            if dimless:
                vel_steps[i + 1, :, :] = vel_steps[i, :, :] + force * timestep
            else:
                vel_steps[i + 1, :, :] = vel_steps[i, :, :] + force * timestep / ARG_MASS

        kin_steps[i+1] = kinetic_energy(vel_steps[i+1, :, :])[1]
        
        energy_steps[0, i+1] = kinetic_energy(vel_steps[i+1, :, :])[1]
        energy_steps[1, i+1] = potential_energy(rel_dis)[2]
        energy_steps[2, i+1] = energy_steps[0, i+1]+energy_steps[1, i+1]

        
        if rescaling and (i < rescaling_max_timesteps+1 and rescale_counter<rescaling_max_rescales and rescaling_moment
                          + average_rescale > i): #add rescaling moment
            # Rescale velocity
            if rescaling_mode == 0:
                v_lambda = 1 # never use this
            else:
                # target kin energy
                target_energy = (num_atoms - 1)*temp*EPSILON*3/2
                current_energy = np.sum(kin_steps[i-average_rescale+1 : i])/average_rescale          #here something nice could have been implemented for efficiency
                need_rescaling = np.abs(target_energy - current_energy) > rescaling_delta_energy

            if need_rescaling:
                rescaling_moment = i
                # limit rescaling factor between 0.5 and 2.0
                v_lambda = np.sqrt(2*target_energy/current_energy)
                #print(v_lambda)
                #print(target_energy)
                #print(current_energy)

                if rescaling_limit:
                    v_lambda = max(rescaling_limit_lower,v_lambda)
                    v_lambda = min(rescaling_limit_upper,v_lambda)

                # apply rescaling factor
                vel_steps[i + 1, :, :] *= v_lambda
                # rescaling statistics below
                rescale_counter+=1
                print(v_lambda)

                need_rescaling = False
                #rescale_max_print = max(rescale_max,v_lambda)
                #rescale_min_print = min(rescale_min,v_lambda)


    #if rescaling:
        # print rescaling statistics
        #print("Rescaled",rescale_counter,"times with λ: [",rescale_min_print,"~",rescale_max_print,"]")
    
    times = np.linspace(0, num_tsteps*dt, num_tsteps)
    plt.plot(times, energy_steps.T)
    plt.xlabel('Time (seconds)')
    plt.ylabel('Energy (dimless)')
    plt.legend(('kinetic energy', 'potential energy', 'total energy'))
    plt.show()

    return pos_steps, vel_steps, energy_steps


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
                    dudt[i, j] = -24 * ((2*(rel_dist[i, j]**-13)) - (rel_dist[i, j]**-7)) / (
                        rel_dist[i, j])
                    
                else:
                    dudt[i, j] = 0
        for i in range(0, len(rel_pos[1])):  # particle i
            for j in range(0, len(rel_pos[1])):  # particle i rel to j (!=i)
                if i!= j:
                    force[i, j, :] = dudt[i, j] * rel_pos[i, j, :]
                else:
                    force[i, j, :] = 0

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
        D: Diffusion coefficient according to lecture notes. NEEDS ELABORATION
        -------

        """

    # make positions continuous
    displacement = 0.0
    # make array with same size
    p00 = np.zeros_like(loc)
    # time iteration
    for k in range(len(loc[0, 0, :])):
        for j in range(len(loc[0, :, 0])):
            for i in range(len(loc[:, 0, 0])):
                p00[i, j, k] = loc[i, j, k] + displacement
                # last value check
                if i != len(loc[:, j, k]) - 1:
                    # check for discontinuity
                    if loc[i + 1, j, k] > loc[i, j, k] + box_dim / 2:
                        displacement -= box_dim
                    if loc[i + 1, j, k] + box_dim / 2 < loc[i, j, k]:
                        displacement += box_dim

    # initiate reference values

    # initiate reference values

    init_loc = p00[timestep, :, :]
    # print('currently we take an average for initial location!
    # see diff coefficient and autocorr function to see if it makes sense')
    # init_loc = np.mean(p00[int(timestep*0.99):int(timestep*1.01), :, :], axis=0)
    print(init_loc)
    loc_usage = p00[timestep:-1, :, :]
    msd_1 = np.abs((loc_usage - init_loc) ** 2)
    # next
    for i in range(len(loc_usage[:, 0, 0])):
        msd_1[i, :, :] = msd_1[i, :, :] / (i + 1)

    msd_2 = np.sum(msd_1, axis=2)
    number_particles = len(program[0][0, :, 0])  # number of particles
    msd_3 = np.sum(msd_2, axis=1) / number_particles
    print(len(msd_3))
    diffusion = np.empty(len(msd_3))
    for i in range(len(msd_3)):
        diffusion[i] = msd_3[i] / (6 * (i + 1))

    plt.plot(p00[:, :, 0])
    plt.title('continous location all particles for x direction')
    plt.show()
    plt.plot(p00[:, :, 1])
    plt.title('continous location all particles for y direction')
    plt.show()
    plt.plot(p00[:, :, 2])
    plt.title('continous location all particles for z direction')
    plt.show()
    plt.plot(msd_2[:, :])
    plt.title('the mean square displacement for each particle summed over all directions')
    plt.show()
    plt.plot(diffusion)
    plt.title('The diffusion coefficient')
    plt.show()
    # print('the diff coeff is shown in the plot above', D)
    return msd_1, msd_2, msd_3, diffusion


# +
# q = ms_displacement(program[0], 15000)
# plt.plot(program[0][:,:,2])
# plt.plot(q[0][:,:,2])
# plt.plot(q[3])

# +
# plt.plot(q[3])
# -

def msd_plot(msd, partnum):
    """"
    plots the MSD of a single atom NOTE MIGHT NOW WORK
    Parameters
    ----------
    msd: np.ndarray
        the msq time dependent array, summed over the dimensions, for M particles
        [msd_part1(dtime=0), msd_part2(dtime=0),... ], [msd_part1(dtime=1), msd_part2(time=1),... ], ....
        best use case: msd_2 from ms_displacement function
    partnum: int
        the particle that is to be plotted by the function
    Returns
    ----------
    None
    """

    plt.plot(msd[:, partnum])
    plt.show()
    return


def auto_corr(data_values, skipvalues):
    """"
    gives the normalized autocorrelation function of an obersvable function.

    Parameters
    ---------------
    data_values: np.ndarray 1D
        The data values used corresponding to the expectation value. This should be an 1D array
        most likely, this is ms_deviation[2]
    skipvalues: int
        skips these initial values. NOTE KEEP AT 0 FOR REPORT.
    Returns
    -------------
    Autocorrelation: np.ndarray
        The autocorrelation function for t
    """
    data_values = data_values[skipvalues:-1]
    number_particles = len(data_values)
    autoc = np.zeros(number_particles)
    for i in range(number_particles - 1):
        nmax = number_particles - i
        ant = data_values[i:number_particles:1]
        an = data_values[0:nmax:1]
        autoc[i] = ((number_particles - i) * np.sum((ant * an)) - (np.sum(an) * np.sum(ant))) / (
                np.sqrt((number_particles - i) * np.sum(an ** 2) - np.sum(an) ** 2) * np.sqrt(
                    (number_particles - i) * np.sum(ant ** 2) - np.sum(ant) ** 2))
    plt.plot(autoc)
    plt.title('The autocorrelation function')
    plt.show()
    return autoc


def exponential_fit(y_data, cutoff):
    """"
    Gives exponential fit of ydata given, removing everything after the cutoff index.
    Note: does not use initial guesses. Check manually from graph if it is okay.

    Parameters
    ---------------
    y_data: np.ndarray 1D
        The data that is to be fitted.
    cutoff: int
        the last datapoint of ydata that is to be used.

    Returns
    -------------
    Params, Tau: float
        fit parameters of the exponential fit
    params_covariance, Covarance of Tau: float
        covariance of tau
    All return values are only taking the data before the cutoff y_data
    """
    numpoints = len(y_data[0:cutoff])
    x_data = np.linspace(0, numpoints, num=numpoints)

    def funcexp(x, tau):
        return np.exp(-x / tau)

    params, params_covariance = optimize.curve_fit(funcexp, x_data, y_data[0:cutoff])
    print('Tau is ', params[0], 'and the covariance of Tau is', params_covariance[0])

    plt.plot(funcexp(x_data, params[0]), 'r', label='fitted autocorrelation')
    plt.plot(y_data[0:cutoff], 'b', label='autocorrelation from data')
    plt.title('The fitted and original autocorrelation function')
    plt.legend()
    plt.show()
    return params, params_covariance


# +
# plt.plot(Q)
# plt.plot(Q[0:4000])
# -
def process_data(positions, velocities):
    print("Test if the total energy is conserved")
    pos1 = positions[0, :, :]
    pos2 = positions[steps - 1, :, :]

    vel1 = velocities[0, :, :]
    vel2 = velocities[steps - 1, :, :]

    r_pos1 = atomic_distances(pos1, box_dim)
    r_pos2 = atomic_distances(pos2, box_dim)

    print("Initial total energy: " + str(total_energy(vel1, r_pos1[1])))
    print("Final total energy:   " + str(total_energy(vel2, r_pos2[1])))
    print("Delta total energy:   " + str(total_energy(vel2, r_pos2[1]) - total_energy(vel1, r_pos1[1])))

    times = np.linspace(0, dt * steps, steps)
    if num_atoms == 2:
        print("Plot inter-atom distance over time")
        distances = np.zeros(steps)
        if dimless:
            for x in range(steps):
                distances[x] = np.max(atomic_distances(positions[x, :, :], box_dim)[1])
            # distances = [np.max(atomic_distances(positions[x, :, :], box_dim)[1]) for x in range(steps)]
        else:
            distances = [np.max(atomic_distances(positions[x, :, :], box_dim)[1]) for x in range(steps)]
        plt.plot(times, distances)
        if dimless:
            plt.ylabel('Distance (dimless)')
            plt.xlabel('Time (dimless)')

        else:
            plt.ylabel('Distance (m)')
            plt.xlabel('Time (s)')

        plt.show()

    print("Print energy levels over time")
    energies = np.zeros([3, steps])
    if dimless:
        for x in range(steps):
            energies[0, x] = kinetic_energy(velocities[x, :, :])[1]
            energies[1, x] = potential_energy(atomic_distances(positions[x, :, :], box_dim)[1])[2]
            energies[2, x] = total_energy(velocities[x, :, :],
                                          atomic_distances(positions[x, :, :], box_dim)[1])
        # energies = [(kinetic_energy(velocities_store[x, :, :])[1],
        #             potential_energy(atomic_distances(positions_store[x, :, :], box_dim)[1])[2],
        #             total_energy(velocities_store[x, :, :],
        #                          atomic_distances(positions_store[x, :, :], box_dim)[1]))
        #           for x in range(steps)]
        # energies = np.array(energies)
    else:
        energies = [kinetic_energy(velocities[x, :, :])[1] for x in range(steps)]

    # times = np.linspace(0, dt*steps, steps)
    plt.plot(times, energies.T)
    plt.xlabel('Time (dimless)')
    plt.ylabel('Energy (dimless)')
    plt.legend(('kinetic energy', 'potential energy', 'total energy'))
    plt.show()
    return energies


def test_initial_velocities(init_velocities):
    if init_velocities is None:
        init_velocities = init_velocity(1000, TEMP, dim)

    vel_mag = np.linalg.norm(init_velocities, axis=1)
    # [np.sqrt(np.sum([v[i] ** 2 for i in range(dim)])) for v in init_velocities]

    plt.hist(vel_mag, bins=15)
    plt.show()


def main():
    """"
        Beginning of program
        to start with "easy" locations and velocities,


    """
    #    easy, handpicked initial positions and velocities.
    # init_pos = [[9.9, 9.6, 8.7], [0.3, 0.6, 0.3], [3.5, 4.6, 5.7], [9.9, 3.3, 6.6], [6.0, 7.5, 9.0],
    #            [0.6, 0.6, 9.0], [3.3, 3.3, 3.3], [8.8, 2.7, 6.3], [6.3, 8.7, 1.5]]
    # init_vel = [[1.2, 0.8, 1.2], [-0.9, -0.67, -0.88], [-0.89, 0.94, 1.55], [1.52, -0.53, 0.97], [0.60, -1.55, 0.32]
    #    , [-0.22, 1.53, -0.34], [1.25, 0.66, -0.97], [-0.36, -1.29, 0.09], [1.22, 0.01, -0.61]]

    # Below is the must be uncommented for the delivarble.
    init_pos = fcc_lattice(num_atoms, box_dim / 2)
    init_vel = init_velocity(num_atoms, TEMP, dim)

    # init_pos = np.array([[1, 0.1, 0.1], [0.1, 0.1, 4.9]])
    # init_vel = np.array([[0.1, 0.1, 0.1], [0.01, 0.01, 0.01]])

    # init_pos = [[25,25,25], [28,25,25], [25,25,27]]
    # init_vel = [[1,0,0], [1,0,0], [1,0,0]]

    test_initial_velocities(None)

    p1, v1, e1 = simulate(init_pos, init_vel, steps, dt, box_dim)
    process_data(p1, v1)

    return p1, v1, e1


program, vel, e1 = main()
print("--- %s seconds ---" % (time.time() - start_time))

# plot the autocorrelation function
q = ms_displacement(program[0], int(rescaling_max_timesteps * 1.2))
focusdiff = 0
plt.title(('The Diffusion coefficient skipping the first', str(focusdiff), 'values'))
plt.plot(q[3][focusdiff:])
plt.show()
qq = auto_corr(q[2], 0)
plotfocus = 500
plt.plot(qq[0:plotfocus])
plt.title(('The autocorrelation function for the first ' + str(plotfocus) + ' values'))
plt.show()
qqq = exponential_fit(qq, plotfocus)
plt.plot(q[1][0:300])
plt.show()

# +
times = np.linspace(0, steps*dt, steps)
plt.plot(times, e1.T)
plt.xlabel('Time (seconds)')
plt.ylabel('Energy (dimless)')
plt.legend(('kinetic energy', 'potential energy', 'total energy'))
plt.show()


# -


