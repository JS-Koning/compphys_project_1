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
            displacement = 0.0

    init_loc = p00[timestep, :, :]
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
    plt.plot(msd_3[:, :])
    plt.title('the mean square displacement')
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


def auto_corr2(data):
    # plot the autocorrelation function (NEEDS TO BE MOVED TO utililities.py)
    q = ms_displacement(data, int(sim.rescaling_max_timesteps * 1.2))
    focusdiff = 0
    plt.title(('The Diffusion coefficient skipping the first', str(focusdiff), 'values'))
    plt.plot(q[3][focusdiff:])
    plt.show()
    qq = auto_corr(q[2], 0)
    plotfocus = 75 #sim.steps/50
    plt.plot(qq[0:plotfocus])
    plt.title(('The autocorrelation function for the first ' + str(plotfocus) + ' values'))
    plt.show()
    exponential_fit(qq, plotfocus)
    plt.plot(q[1][0:300])
    plt.show()


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

    print("Initial total energy: " + str(sim.total_energy(vel1, r_pos1[1])[0]))
    print("Final total energy:   " + str(sim.total_energy(vel2, r_pos2[1])[0]))
    print("Delta total energy:   " + str(sim.total_energy(vel2, r_pos2[1])[0] - sim.total_energy(vel1, r_pos1[1])[0]))

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
            t, k, p = sim.total_energy(velocities[x, :, :],
                                          sim.atomic_distances(positions[x, :, :], sim.box_dim)[1])
            energies[0, x] = k
            energies[1, x] = p
            energies[2, x] = t
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


def test_initial_velocities(init_velocities):
    if init_velocities is None:
        init_velocities = sim.init_velocity(1000, sim.TEMP, sim.dim)

    vel_mag = np.linalg.norm(init_velocities, axis=1)
    # [np.sqrt(np.sum([v[i] ** 2 for i in range(dim)])) for v in init_velocities]

    gaussian_mean = np.mean(vel_mag)
    gaussian_sigma = np.std(vel_mag)**2
    gaussian_max = np.max(vel_mag)

    x_axis = np.linspace(0.5,3.0,1000)

    gaussian = np.exp(-np.power(x_axis-gaussian_mean,2.0)/gaussian_sigma/2)

    y, x, _ = plt.hist(vel_mag, bins=15)

    gaussian *= np.max(y)
    plt.plot(x_axis, gaussian)
    plt.show()


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
    data_values = data_values[skipvalues:]
    N = len(data_values)
    autoc = np.zeros(N)
    for i in range(N-1):
        nmax = N - i
        Ant = data_values[i:N:1]
        An = data_values[0:nmax:1]
        autoc [i] = ((N-i) * np.sum( (Ant*An) ) - ( np.sum(An) * np.sum(Ant) )) / ( np.sqrt((N-i) * np.sum(An**2) - np.sum(An)**2) * np.sqrt((N-i) * np.sum(Ant**2) - np.sum(Ant)**2)  )
    plt.plot(autoc)
    plt.title('The autocorrelation function')
    plt.show()
    return autoc


def exponential_fit(y_data, cutoff):
    """"
    Gives exponential fit of ydata given, removing everything after the cutoff index. Note: does not use initial guesses. Check manually from graph if it is okay.
    
    Parameters
    ---------------
    ydata: np.ndarray 1D
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
    x_data = np.linspace(0,numpoints, num=numpoints)
    def funcexp(x, tau):
        return np.exp(-x/tau)
    
    params, params_covariance = optimize.curve_fit(funcexp, x_data, y_data[0:cutoff])
    print('Tau is ', params[0], 'and the covariance of Tau is', params_covariance[0])
    
    plt.plot(funcexp(x_data, params[0]), 'r', label='fitted autocorrelation')
    plt.plot(y_data[0:cutoff], 'b', label='autocorrelation from data')
    plt.title('The fitted and original autocorrelation function')
    plt.legend()
    plt.hlines(0, 0, cutoff)
    plt.show()
    return params, params_covariance


def expectedvalues(y_data, cutoff):
    """"
    Gives expected value according to <A> = 1/N * sum(An) with n>0 as lower boundary, and N as upper boundary.
    Please note, here the index starts at 0, not at 1 as in the literature.
    
    Parameters
    ---------------
    y_data: np.ndarray 1D
        The input array, An.
    cutoff: int
        the last datapoint of y_data that is to be used.
    
    Returns
    -------------
    expected: float
        The expected value of y_data <A>
    squared_expected: float
        The squared expected value <A**2>
    expected_squared
        the squared value of the expected value <A>**2    
    """
    A = y_data[:cutoff]
    N = len(A)

    # Python 3.9 fix... (works fine in Python 3.8)
    #for x in range(len(A)):
    #    val = A[x]
    #    if (np.abs(val) > 1) or val == float("inf") or val != val:
    #        print("Error:", val)
    #        A[x] = 0

    sumA = np.sum(A)

    expected = 1/N * sumA
    expected2 = np.power(expected,2.0)
    square_expected = 1/N * sum(A**2)

    return expected, expected2, square_expected


def errortau(y_data, tau):
    """"
    calculates the error in the mean of the autocorrelation function.
    
    Parameters
    ---------------
    ydata: np.ndarray 1D
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
    N = len(y_data)
    expectedA = expectedvalues(y_data, N)
    sigma = expectedA[2] - expectedA[1]
    sigmaA = np.sqrt(2*tau*sigma/N)
    return sigmaA, sigma


def block_data(y_data, block_length):
    """"
    Takes average of the block_length values, and puts this in a array.
    
    Parameters
    ---------------
    y_data: np.ndarray 1D
        The input array that requires data blocking.
    block_length: integer
        The required block length.
    
    Returns
    -------------
    a: np.ndarray 1D
        The new array of the block, created by taking the block averaged of the y_data.
    """
    Nb = len(y_data)//block_length
    a = np.empty(Nb)
    for i in range(1,Nb):
        a[i] = sum(y_data[((i-1)*block_length)+1:i*block_length])/block_length
    np.delete(a, 0)
    return a


def errorblock(meanblocks):
    """"
    Calculated the error of the mean, taking the datablocks as input. These datablocks can be found using the block_data function.
    Note: This is dependent on the size of the block!
    
    Parameters
    ---------------
    meanblocks: np.ndarray 1D
        Input array, in the literature of this course called a_i
    y_data: np.ndarray 1D
        The input array that requires data blocking.
    block_length: integer
        The required block length.
    
    Returns
    -------------
    sigmaAb: np.ndarray 1D
        the standard deviation of the estimator of the mean (error of the mean).
    """
    expecteda = expectedvalues(meanblocks, len(meanblocks))
    sigmaAb = np.sqrt((expecteda[2]-expecteda[1])/(len(meanblocks)-1))
    return sigmaAb


def error_mean(y_data, cutoff):
    """"
    Calculates the error in the mean of the observable
    Shows the plot of the autocorrelation function.
    To verify results:
    Check where the errorvsblocksize converges.
    
    Parameters
    ---------------
    y_data: np.array 1D
        observable
    cutoff: int
        the cutoff value determined by the autocorrelation function.
    
    Returns
    -------------
    none
    """
    
    y_data = ms_displacement(y_data, 0)[3]
    autofun = auto_corr(y_data, 0)
    fit = exponential_fit(autofun, cutoff)
    max_block_size = int(len(y_data)/20)
    errora = np.empty(max_block_size)
    for i in range(2,max_block_size):
        blocks = block_data(y_data, i)
        errora[i] = errorblock(blocks)
    plt.plot(errora)
    plt.title('error vs block size')
    plt.show()
    tauer = errortau(y_data, fit[0])
    print('uncertainty in the mean is of the mean squared distance ', tauer[0])
    return



