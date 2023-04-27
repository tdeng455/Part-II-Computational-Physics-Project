"""initial.py
Creates a lattice and defines useful functions such as calculating neighbouring sites, 
calculating the sum of neighbouring spin sites and computes the total magnetisation of the lattice
"""

import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
from statsmodels.tsa.stattools import acf
import types

def create_lattice(width,type=0):
    """Creates a lattice array of three types
    
    Type 1: all spins +1
    Type -1: all spins -1
    Type 0: a random lattice of spins +/-
    Type 2: anti-aligned lattice
    """
    lattice = np.array((width,width))   
    if type == 1:
        return np.ones((width,width))
    elif type == 0:
        return np.random.choice([1,-1],(width,width))
    elif type == -1:
        return (-1)*np.ones((width,width))
    elif type == 2:
            if width%2 == 0:
                return np.tile([[1,-1],[-1,1]],(width//2,width//2))
            else:
                return np.tile([[1,-1],[-1,1]],((width+1)//2,(width+1)//2))[:width,:width]
    else:
        raise ValueError("Invalid type. Type should be 0, 1, or -1.")

def plot_lattice(lattice,ax,title):
    """Plots the lattice with black as +1, and white as -1"""
    ax.matshow(lattice, vmin=-1, vmax=1, cmap=plt.cm.binary)
    ax.title.set_text(title)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_xticks([])

def get_neighbouring_sites(i,j,width):
    """Returns the coordinates of the neighbouring sites of a given coordinate with periodic boundrary conditions"""
    neighbours = [((i-1)%width, j), ((i+1)%width, j), (i, (j-1)%width), (i, (j+1)%width)]
    return neighbours

def neighbouring_spins_sum(i,j,lattice,width):
    """Returns the sum of the spins of the neighbouring sites of a given coordinate (i,j)"""
    neighbours = get_neighbouring_sites(i,j,width)
    sum_neighbours = sum(lattice[i] for i in neighbours)
    return sum_neighbours

def magnetisation(lattice):
    """Computes the overall magnetisation of the lattice"""
    Mag = lattice.sum()/np.size(lattice)
    return Mag

"""Defining autocorrelation function and autocorrelation times 
   using statsmodels library
"""
def autocorrelation(data):
    N = len(data)
    autocorr = acf(data,adjusted=False,fft=False, nlags=(N-1))
    return autocorr

def autocorrelation_time(data):
    autocorr = autocorrelation(data)
    crit = np.exp(-1)
    time = np.argmin(autocorr>crit, axis=0)
    return time

def batch_estimate(data, operation, num_batches, batch_with_autocorr):
    if batch_with_autocorr == True:
        t_f = autocorrelation_time(data)
        print(t_f)
        if t_f == 0:
            m = int(len(data)/2)
        else: 
            m = int(len(data)/(2*t_f))

    elif batch_with_autocorr == False:
        m = num_batches

    else:
        raise ValueError("Invalid entry. batch_with_autocorr should be True or False")
    
    if isinstance(operation, types.FunctionType) == False:
        raise ValueError("Invalid entry, operation variable must be a function on an array")
    
    batches = np.array_split(data, m)
    assign = np.array(list(map(operation, batches)))
    estimate = np.mean(assign)
    error = np.std(assign)/np.sqrt(m-1)
    return estimate, error

#############################################################################

"""Defining autocorrelation function and autocorrelation times
   without statsmodels library
"""

def autocorrelation_2(data,t_max):
    """Computes ACF for a given time series"""
    size = len(data)
    mean = np.mean(data)

    # autocovariance
    autocov = np.zeros(size)
    for t in range(min(t_max,size)):
        autocov[t] = np.dot(data[:size-t] - mean, data[t:] - mean) / (size)

    #normalise
    autocorr = autocov/autocov[0]

    return autocorr

def autocorrelation_time_2(data, tmax):
    """Finds the autocorrelation time of a given data series"""
    autocorr = autocorrelation_2(data, tmax)
    crit = np.exp(-1)
    time = np.argmin(autocorr>crit, axis=0)
    return time