"""Calculates and plots autocorrelation function using statsmodels library"""

import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
import functions.initial as initial
import functions.metropolis as metropolis
import functions.cluster as cluster

width = 40  # width of the lattice
betaJ = 0.01

# Metropolis-Hastings ACF
print('MH mags')
mag_MH = []
MH_n_steps = 500  # number of steps to run each simulation for Metropolis-Hastings algorithm
MH_burn_in = 100000  # euilibration time for Metropolis-Hastings algorithm

lattice = initial.create_lattice(width, 0)
metropolis.n_MH_moves(lattice, width, betaJ, MH_burn_in) # equilibration

for step in range(MH_n_steps):
    metropolis.MH_flip(lattice, width, betaJ)
    mag_MH.append(abs(initial.magnetisation(lattice)))
    if step%10 == 0:
        print(step)

#Wolf ACF
print('Wolff mags')
p_add = 1 - np.exp(-2*betaJ)
mag_wolff = []
wolff_n_steps = 500 # number of moves for wolff algorithm
wolff_burn_in = 100 # equilibration time for wolff

lattice = initial.create_lattice(width, 0)
cluster.n_wolff_moves(lattice, p_add, wolff_burn_in) # equilibration

for step in range(wolff_n_steps):
    cluster.wolff_flip1(lattice, p_add)
    mag_wolff.append(abs(initial.magnetisation(lattice)))
    if step%10 == 0:
        print(step)

print('Calc of Metropolis-Hastings ACF')
ACF_MH = initial.autocorrelation(mag_MH)
print('Calc of Wolff ACF')
ACF_wolff = initial.autocorrelation(mag_wolff)

#data = [ACF_MH,ACF_wolff]
#np.save('acf_samemoves4', data)


plt.figure()
plt.plot(ACF_wolff, label='Wolff', color='blue')
plt.plot(ACF_MH, label='Metropolis-Hastings', color='red')
plt.axhline(y=np.exp(-1), color='black', linestyle='--', label='1/e')
plt.xlabel('Lag time')
plt.ylabel('ACF')
plt.title('ACF for MH vs Wolff at beta = {}'.format(betaJ))
plt.legend()
plt.show()