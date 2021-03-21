import simulate as sim
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

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
                    if loc[i + 1, j, k] > loc[i, j, k] + sim.box_dim / 2:
                        displacement -= sim.box_dim
                    if loc[i + 1, j, k] + sim.box_dim / 2 < loc[i, j, k]:
                        displacement += sim.box_dim

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
    number_particles = len(loc[0, :, 0])  # number of particles
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
    pos2 = positions[sim.steps - 1, :, :]

    vel1 = velocities[0, :, :]
    vel2 = velocities[sim.steps - 1, :, :]

    r_pos1 = sim.atomic_distances(pos1, sim.box_dim)
    r_pos2 = sim.atomic_distances(pos2, sim.box_dim)

    print("Initial total energy: " + str(sim.total_energy(vel1, r_pos1[1])))
    print("Final total energy:   " + str(sim.total_energy(vel2, r_pos2[1])))
    print("Delta total energy:   " + str(sim.total_energy(vel2, r_pos2[1]) - sim.total_energy(vel1, r_pos1[1])))

    times = np.linspace(0, sim.dt * sim.steps, sim.steps)
    if sim.num_atoms == 2:
        print("Plot inter-atom distance over time")
        distances = np.zeros(sim.steps)
        if sim.dimless:
            for x in range(sim.steps):
                distances[x] = np.max(sim.atomic_distances(positions[x, :, :], sim.box_dim)[1])
            # distances = [np.max(atomic_distances(positions[x, :, :], box_dim)[1]) for x in range(steps)]
        else:
            distances = [np.max(sim.atomic_distances(positions[x, :, :], sim.box_dim)[1]) for x in range(sim.steps)]
        plt.plot(times, distances)
        if sim.dimless:
            plt.ylabel('Distance (dimless)')
            plt.xlabel('Time (dimless)')

        else:
            plt.ylabel('Distance (m)')
            plt.xlabel('Time (s)')

        plt.show()

    print("Print energy levels over time")
    energies = np.zeros([3, sim.steps])
    if sim.dimless:
        for x in range(sim.steps):
            energies[0, x] = sim.kinetic_energy(velocities[x, :, :])[1]
            energies[1, x] = sim.potential_energy(sim.atomic_distances(positions[x, :, :], sim.box_dim)[1])[2]
            energies[2, x] = sim.total_energy(velocities[x, :, :],
                                          sim.atomic_distances(positions[x, :, :], sim.box_dim)[1])
        # energies = [(kinetic_energy(velocities_store[x, :, :])[1],
        #             potential_energy(atomic_distances(positions_store[x, :, :], box_dim)[1])[2],
        #             total_energy(velocities_store[x, :, :],
        #                          atomic_distances(positions_store[x, :, :], box_dim)[1]))
        #           for x in range(steps)]
        # energies = np.array(energies)
    else:
        energies = [sim.kinetic_energy(velocities[x, :, :])[1] for x in range(sim.steps)]

    # times = np.linspace(0, dt*steps, steps)
    plt.plot(times, energies.T)
    plt.xlabel('Time (dimless)')
    plt.ylabel('Energy (dimless)')
    plt.legend(('kinetic energy', 'potential energy', 'total energy'))
    plt.show()
    return energies


def test_initial_velocities(init_velocities):
    if init_velocities is None:
        init_velocities = sim.init_velocity(1000, sim.TEMP, sim.dim)

    vel_mag = np.linalg.norm(init_velocities, axis=1)
    # [np.sqrt(np.sum([v[i] ** 2 for i in range(dim)])) for v in init_velocities]

    plt.hist(vel_mag, bins=15)
    plt.show()