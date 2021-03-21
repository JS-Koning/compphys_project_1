import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import time
# LOAD SEED
from pickle import load

import simulate as sim
import utilities as utils

with open('state.obj', 'rb') as f:
    np.random.set_state(load(f))


def main():
    """"
        Beginning of program
    """

    # Start timer for program run-time
    start_time = time.time()

    # Compare Verlet vs. Euler using simple simulation
    # Simulation parameters
    sim.dim = 3
    sim.box_dim = 3
    sim.num_atoms = 2
    sim.rescaling = False

    # Simple positions and velocities
    init_pos = np.array([[+0.30, +0.45, +0.60], [+2.70, +2.70, +2.70]])
    init_vel = np.array([[-0.15, -0.30, -0.45], [+0.30, +0.45, +0.60]])

    # Verlet
    sim.verlet = True
    pv, vv = sim.simulate(init_pos, init_vel, sim.steps, sim.dt, sim.box_dim)
    utils.process_data(pv, vv)

    # Euler
    sim.verlet = False
    pe, ve = sim.simulate(init_pos, init_vel, sim.steps, sim.dt, sim.box_dim)
    utils.process_data(pe, ve)

    # Keep Verlet for further simulations
    sim.verlet = True
    # Simulation parameters
    sim.dim = 3
    sim.box_dim = 3.41995
    sim.num_atoms = 32
    sim.rescaling = True
    sim.temp = 0.5 * sim.EPSILON / sim.KB

    # For 32 particles
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

    # FCC lattice positions and random velocities
    init_pos = sim.fcc_lattice(sim.num_atoms, sim.box_dim / 2)
    init_vel = sim.init_velocity(sim.num_atoms, sim.TEMP, sim.dim)

    # Test if random velocities are in compliance with Maxwell-Boltzmann distribution
    utils.test_initial_velocities(init_vel)

    #
    p1, v1 = sim.simulate(init_pos, init_vel, sim.steps, sim.dt, sim.box_dim)
    utils.process_data(p1, v1)

    utils.auto_corr2(p1)

    # do tests here

    # End timer for program run-time
    print("--- %s seconds ---" % (time.time() - start_time))
    return p1, v1


main()
