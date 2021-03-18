import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


def normal_autocorr(mu, sigma, tau, N):
    """Generates an autocorrelated sequence of Gaussian random numbers.
    
    Each of the random numbers in the sequence of length `N` is distributed
    according to a Gaussian with mean `mu` and standard deviation `sigma` (just
    as in `numpy.random.normal`, with `loc=mu` and `scale=sigma`). Subsequent
    random numbers are correlated such that the autocorrelation function
    is on average `exp(-n/tau)` where `n` is the distance between random
    numbers in the sequence.
    
    This function implements the algorithm described in
    https://www.cmu.edu/biolphys/deserno/pdf/corr_gaussian_random.pdf
    
    Parameters
    ----------
    
    mu: float
        mean of each Gaussian random number
    sigma: float
        standard deviation of each Gaussian random number
    tau: float
        autocorrelation time
    N: int
        number of desired random numbers
    
    Returns:
    --------
    sequence: numpy array
        array of autocorrelated random numbers
    """
    f = np.exp(-1./tau)
    
    sequence = np.zeros(shape=(N,))
    
    sequence[0] = np.random.normal(0, 1)
    for i in range(1, N):
        sequence[i] = f * sequence[i-1] + np.sqrt(1 - f**2) * np.random.normal(0, 1)
    
    return mu + sigma * sequence


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


def errorauto(y_data, tau):
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
    N = len(y_data)
    squared = (1/N) * sum((y_data)**2)
    mean = (1/N) * sum(y_data)
    mean2 = mean**2
    var = squared-mean
    error2 = (2*tau/N) * (var)
    sigmaA = np.sqrt(error2)
    print('sigmaA is',sigmaA, 'var is', var, 'mean is ', mean)
    return(sigmaA, var)


def expectedvalues(y_data, cutoff):
     """"
    Gives expected value according to <A> = 1/N * sum(An) with n>0 as lower boundary, and N as upper boundary.
    
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
    A = y_data[1:cutoff]
    N = len(A)
    expected = 1/N * sum(A)
    squared_expected = expected**2
    expected_squared = 1/N * sum(A**2)
    return(expected, squared_expected, expected_squared)


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
    expected: float
        The expected value of y_data <A>
    squared_expected: float
        The squared expected value <A**2>
    expected_squared
        the squared value of the expected value <A>**2    
    """
    Nb = len(y_data)//block_length
    a = np.empty(Nb)
    for i in range(1, Nb):
        a[i] = sum(y_data[(i-1)*block_length+1:i*block_length])


def block_data(y_data, tau)
    tau = int(tau)
    Nb = int(len(y_data)/(2*tau))
    blocksize = 2*tau
    blockavg = np.empty(Nb)
    for i in range(1,Nb):
        blockavg[i] = np.average(y_data[i*blocksize+1:(i+1)*blocksize:1])
    return(blockavg)


#data block error(b)
def errora(block_data):
    block_data = block_data[1:]
    Nb = len(block_data)
    squared = (1/Nb) * sum(block_data**2)
    mean = (1/Nb) * sum(block_data)
    mean2 = mean**2
    error2 = (1/(Nb-1)) * (squared-mean2)
    #error2 = squared-mean2
    error = np.sqrt(error2)
    return(error)


# +
tau = 50
mu = 0
sigma = 1
N = 20000


y_data = normal_autocorr(mu, sigma, tau, N)
auto = auto_corr(y_data, 0)
plt.plot(y_data)
plt.show()


# -

autocorr = auto_corr(y_data, 0)
a = exponential_fit(autocorr, 150)

taufound = 41.8476123563759
skimmed = y_data[0:int(taufound)]
leng = len(skimmed)
print(skimmed)
print(sum(skimmed))
a = sum(skimmed)
b = sum(skimmed/leng)
print(b)
skimmed[0]

xax = 290
er = np.empty(xax)
for i in range(1,xax):
    q = block_data(y_data, i)
    er[i] = errora(q)
plt.plot(er)


errorauto(y_data, 41.8476123563759)


