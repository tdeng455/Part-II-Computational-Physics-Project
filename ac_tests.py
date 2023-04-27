"""RUNTIME TESTS FOR AC FUNCTIONS"""

import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
import functions.initial as initial
import functions.metropolis as metropolis
import functions.cluster as cluster
import time

temps = np.arange(1.5,3.5,0.5)  # temperatures to simulate
width = 40  # width of the lattice


# Metropolis-Hastings algorithm simulation
mags_MH = []
MH_n_steps = 150000  # number of steps to run each simulation for Metropolis-Hastings algorithm
MH_burn_in = 50000  # equilibration time for Metropolis-Hastings algorithm
for temp in temps:
    print('temp = ' + str(temp))
    betaJ = 1/temp
    lattice = initial.create_lattice(width, 1)
    ms = []
    metropolis.n_MH_moves(lattice, width, betaJ, MH_burn_in)
    
    for step in range(MH_n_steps):
    
        metropolis.n_MH_moves(lattice, width, betaJ, 1)
    
        if step%10 == 0:
            m = np.abs(initial.magnetisation(lattice))
            ms.append(m)
    
    mags_MH.append(ms)

"""
# Wolff algorithm simulation
mags_wolff = []
wolff_n_steps = 1000  # number of steps to run each simulation for Wolff algorithm
wolff_burn_in = 100  # equlibration time for Wolff algorithm
for temp in temps:
    print('temp = ' + str(temp))
    p_add = 1-np.exp(-2/temp)
    lattice = initial.create_lattice(width, 1)
    ms = []
    cluster.n_wolff_moves(lattice, p_add, wolff_burn_in)
    
    for step in range(wolff_n_steps):
        cluster.n_wolff_moves(lattice, p_add, 1)
        #print(step)
        m = np.abs(initial.magnetisation(lattice))
        ms.append(m)
    
    mags_wolff.append(ms)
"""

#Plot and save autocorrelation time against temperature

MH_autocorr_times = []
for mags in mags_MH:
    print('(slower?) start')
    start_time = time.time()
    MH_autocorr_times.append(initial.autocorrelation_time_2(mags,len(mags)))
    end_time = time.time()
    print('(slower?) end, t = ' + str(end_time - start_time))

MH_autocorr_times_fast = []
for mags in mags_MH:
    print('(faster?) start')
    start_time1 = time.time()
    MH_autocorr_times_fast.append(initial.autocorrelation_time(mags))
    end_time1 = time.time()
    print('(faster?) end, t = ' + str(end_time1 - start_time1))

#save
#MH_autocorr_dat = [temps,MH_autocorr_times]
#np.save('MH_autocorrelation_data', MH_autocorr_dat)

#wolff_autocorr_times = []
#for mags in mags_wolff:
#    print(len(mags))
#    wolff_autocorr_times.append(initial.autocorrelation_time(mags))

#save
#wolff_autocorr_dat = [temps,wolff_autocorr_times]
#np.save('wolff_autocorrelation_data', wolff_autocorr_dat)

fig, ax = plt.subplots(1,2)
ax[0].plot(temps, MH_autocorr_times, label='Metropolis-Hastings slow')
ax[1].plot(temps, MH_autocorr_times_fast, label='Metropolis-Hastings fast')
#plt.plot(temps, wolff_autocorr_times, label='Wolff')
ax[0].set_xlabel('Temperature')
ax[0].set_ylabel('Autocorrelation Time')
ax[1].set_xlabel('Temperature')
ax[1].set_ylabel('Autocorrelation Time')
fig.suptitle('Autocorrelation Time vs Temperature for each Algorithm')
plt.legend()
plt.show()