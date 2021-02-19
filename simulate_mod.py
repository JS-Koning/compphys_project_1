from constants_m import *
import numpy as np

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
    atomic_distances_step = np.zeros(num_tsteps, num_atoms, dim)

    pos_steps[0, :, :] = init_pos
    vel_steps[0, :, :] = init_vel
    for i in range(num_tsteps-1):
        # Function that applies periodic bc IFF required
        # Note! this is not yet written explicitly
        pos_steps[i, :, :] = periodic_func(pos_steps)
        # calculates next position at i+1, particle will  get corrected
        # for box dimensions in next iteration
        pos_steps[i+1, :, :] = pos_steps[i, :, :] + vel_steps[i, :, :]/mass
        # Next we calculate relative distance of these particles
        atomic_distances_step[i, :, :] = atomic_distances(pos_steps[i, :, :], box_dim)
        # calculate next velocity, using a undefined, new function that is missing.
        vel_steps[i+1, :, :] = vel_steps[i, :, :] + lj_force(parameters_min_image_pos,parameters_min_image_dist) / mass
        # Potential energy also depends on image method [time, total energy]
        pot_energy[i, :] = potential_energy(parameters_min_image)
        # Kinetic energy does not depend on image method [time, total energy]
        kin_energy[i, :] = kinetic_energy(vel_steps[i, :, :])
    global positions_store
    positions_store = pos_steps
    global velocities_store
    velocities_store = vel_steps

    return pos_steps, vel_steps
