import numpy as np
from constants_m import *
from init_position_m import init_position
from kinetic_energy_m import kinetic_energy
from lj_force_m import lj_force
from atomic_distances import atomic_distances
from potential_energy_m import potential_energy
from simulate_mod import simulate


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

