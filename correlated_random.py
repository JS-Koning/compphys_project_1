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
    for x in range(len(A)):
        val = A[x]
        if (np.abs(val) > 1) or val == float("inf") or val != val:
            print("Error:", val)
            A[x] = 0

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


# +
tau = 50
mu = 0
sigma = 1
N = 20000


y_data = normal_autocorr(mu, sigma, tau, N)
autofun = auto_corr(y_data, 0)
plt.plot(y_data)
plt.show()


# -

fit = exponential_fit(autofun, 300)

max_block_size=300
errora = np.empty(max_block_size)
for i in range(2,max_block_size):
    blocks = block_data(y_data, i)
    errora[i] = errorblock(blocks)
plt.plot(errora)
plt.show()
tauer = errortau(y_data, fit[0])
print(tauer)


