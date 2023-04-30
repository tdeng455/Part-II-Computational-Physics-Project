import numpy as np
import matplotlib.pyplot as plt

# Define the function that generates the Ising model configuration
def initialstate(N):
    ''' Generates a random spin configuration for initial condition'''
    state = 2*np.random.randint(2, size=(N,N))-1
    return state

# Define the function that computes the energy of a spin configuration
def energy(state):
    '''Energy of a given configuration'''
    N = len(state)
    energy = 0
    for i in range(N):
        for j in range(N):
            S = state[i,j]
            nb = state[(i+1)%N, j] + state[i,(j+1)%N] + state[(i-1)%N, j] + state[i,(j-1)%N]
            energy += -nb*S
    return energy/4.

# Define the function that performs one Monte Carlo step
def mcmove(config, beta):
    '''Monte Carlo move using Metropolis algorithm '''
    N = len(config)
    for i in range(N):
        for j in range(N):
            a = np.random.randint(0, N)
            b = np.random.randint(0, N)
            s = config[a, b]
            nb = config[(a+1)%N, b] + config[a,(b+1)%N] + config[(a-1)%N, b] + config[a,(b-1)%N]
            cost = 2*s*nb
            if cost < 0:
                s *= -1
            elif np.random.rand() < np.exp(-cost*beta):
                s *= -1
            config[a, b] = s
    return config

# Define the function that performs Monte Carlo simulation
def simulate(beta, N, num_iter):
    '''Perform Monte Carlo simulation using Metropolis algorithm'''
    config = initialstate(N)
    E = energy(config)
    E_avg = 0
    M_avg = 0
    for i in range(num_iter):
        config = mcmove(config, beta)
        E = energy(config)
        M = np.sum(config)
        E_avg += E
        M_avg += M
    E_avg /= num_iter
    M_avg /= num_iter
    return E_avg, M_avg

# Define the function that performs finite size scaling analysis
def finite_size_scaling(beta_min, beta_max, num_betas, N_min, N_max, num_Ns, num_iter):
    '''Perform finite size scaling analysis to determine the dynamic exponent'''
    beta_vals = np.linspace(beta_min, beta_max, num_betas)
    N_vals = np.linspace(N_min, N_max, num_Ns, dtype=int)
    E_avg = np.zeros((num_betas, num_Ns))
    M_avg = np.zeros((num_betas, num_Ns))
    for j in range(num_Ns):
        N = N_vals[j]
        for i in range(num_betas):
            beta = beta_vals[i]
            E_avg[i,j], M_avg[i,j] = simulate(beta, N, num_iter)
    E_var = np.var(E_avg, axis=1)
    M_var = np.var(M_avg, axis=1)
    E_max = np.argmax(E_var)
    M_max = np.argmax(M_var)
    z = (np.log(N_vals[E_max])-np.log(N_vals[M_max]))/(np.log(beta_vals[E_max])-np.log(beta_vals[M_max]))
    return z

beta_min = 0.1
beta_max = 5
num_betas = 50
N_min = 4
N_max = 64
num_Ns = 5
num_iter = 10000

z = finite_size_scaling(beta_min, beta_max, num_betas, N_min, N_max, num_Ns, num_iter)