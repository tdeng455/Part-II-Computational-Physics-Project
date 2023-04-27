"""Calculates and plots autocorrelation function using acfast"""

import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
import functions.initial as initial
import functions.metropolis as metropolis
import functions.cluster as cluster

width = 40  # width of the lattice

# Metropolis-Hastings ACF
mag_MH = []
MH_n_steps = 20000  # number of steps to run each simulation for Metropolis-Hastings algorithm
MH_burn_in = 10000  # euilibration time for Metropolis-Hastings algorithm
betaJ = 0.4407

lattice = initial.create_lattice(width, 0)
metropolis.n_MH_moves(lattice, width, betaJ, MH_burn_in) # equilibration

for step in range(MH_n_steps):
    metropolis.MH_flip(lattice, width, betaJ)
    mag_MH.append(initial.magnetisation(lattice))
    if step%100 == 0:
        print(step)

#Wolf ACF
p_add = 1 - np.exp(-2*betaJ)
mag_wolff = []
wolff_n_steps = 20000 # number of moves for wolff algorithm
wolff_burn_in = 10000 # equilibration time for wolff

lattice = initial.create_lattice(width, 0)
cluster.n_wolff_moves(lattice, p_add, wolff_burn_in) # equilibration

for step in range(wolff_n_steps):
    cluster.wolff_flip1(lattice, p_add)
    mag_wolff.append(initial.magnetisation(lattice))
    if step%100 == 0:
        print(step)

ACF_MH = initial.autocorrelation(mag_MH)
ACF_wolff = initial.autocorrelation(mag_wolff)

#data = [ACF_MH,ACF_wolff]
#np.save('acf_samemoves2', data)


plt.figure()
plt.plot(ACF_wolff, label='Wolff', color='blue')
plt.plot(ACF_MH, label='Metropolis-Hastings', color='red')
plt.xlabel('Lag time')
plt.ylabel('ACF')
plt.title(r'ACF for MH vs Wolff at critical temperature')
plt.legend()
plt.show()