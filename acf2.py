"""Calculates and plots autocorrelation function against algorithm moves"""

import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
import functions.initial as initial
import functions.metropolis as metropolis
import functions.cluster as cluster

width = 40  # width of the lattice
betaJ = 0.4407

# Metropolis-Hastings ACF
mag_MH = []
MH_n_steps = 2000  # number of steps to run each simulation for MH
MH_burn_in = 10000  # equilibration time for MH

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
wolff_n_steps = 2000 # number of moves for wolff algorithm
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
ACF_MH_2 = initial.autocorrelation_2(mag_MH,len(mag_MH))
ACF_wolff_2 = initial.autocorrelation_2(mag_wolff, len(mag_wolff))

data = [ACF_MH,ACF_wolff]
np.save('acf_comp2', data)
data2 = [ACF_MH_2,ACF_wolff_2]
np.save('acf_basic2', data2)


plt.figure()
plt.plot(ACF_wolff, label='Wolff', color='blue')
plt.plot(ACF_MH, label='Metropolis-Hastings', color='red')
plt.plot(ACF_wolff_2, label='Wolff2', color='blue')
plt.plot(ACF_MH_2, label='Metropolis-Hastings2', color='red')
plt.xlabel('Lag time')
plt.ylabel('ACF')
plt.title(r'ACF for MH vs Wolff at critical temperature')
plt.legend()
plt.show()