# initalizing self defined system parameters
num_atoms = 2  # amount of particles
dim = 3  # dimensions
box_dim = 10  #2* 1.547 #10  # meters; bounding box dimension
dt = 1e-4  # s; stepsize
steps = 30000  # amount of steps
dimless = True  # use dimensionless units
periodic = True  # use periodicity
verlet = True  # use Verlet's algorithm (false: Euler's algorithm)
rescaling = False # use Temperature rescaling
rescaling_mode = 1 # 0 = kin-NRG-based | temp-based
rescaling_delta = 0.0027 # delta for activation of rescaling
rescaling_timesteps = steps / 30 # timesteps interval for rescaling check
rescaling_max_timesteps = steps/2 # max timesteps for rescaling