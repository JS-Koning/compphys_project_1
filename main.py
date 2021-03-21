import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import time
# LOAD SEED
from pickle import load

import simulate as sim
import utilities as utils

start_time = time.time()

with open('state.obj', 'rb') as f:
    np.random.set_state(load(f))

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

    # Below is the random data (must be uncommented for the deliverable)
    init_pos = sim.fcc_lattice(sim.num_atoms, sim.box_dim / 2)
    init_vel = sim.init_velocity(sim.num_atoms, sim.TEMP, sim.dim)

    # init_pos = np.array([[1, 0.1, 0.1], [0.1, 0.1, 4.9]])
    # init_vel = np.array([[0.1, 0.1, 0.1], [0.01, 0.01, 0.01]])

    # init_pos = [[25,25,25], [28,25,25], [25,25,27]]
    # init_vel = [[1,0,0], [1,0,0], [1,0,0]]

    utils.test_initial_velocities(None)

    p1, v1 = sim.simulate(init_pos, init_vel, sim.steps, sim.dt, sim.box_dim)
    utils.process_data(p1, v1)

    return p1, v1


program = main()
print("--- %s seconds ---" % (time.time() - start_time))

# plot the autocorrelation function (NEEDS TO BE MOVED TO utililities.py)
q = utils.ms_displacement(program[0], int(sim.rescaling_max_timesteps * 1.2))
focusdiff = 0
plt.title(('The Diffusion coefficient skipping the first', str(focusdiff), 'values'))
plt.plot(q[3][focusdiff:])
plt.show()
qq = utils.auto_corr(q[2], 0)
plotfocus = 500
plt.plot(qq[0:plotfocus])
plt.title(('The autocorrelation function for the first ' + str(plotfocus) + ' values'))
plt.show()
qqq = utils.exponential_fit(qq, plotfocus)
plt.plot(q[1][0:300])
plt.show()
