import numpy as np
import matplotlib.pyplot as plt
from init_position_m import init_position
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


print(init_position(2,2,2))


