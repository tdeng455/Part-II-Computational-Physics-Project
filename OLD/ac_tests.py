"""RUNTIME TESTS FOR AC FUNCTIONS"""

"""
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
"""

"""Calculates and plots autocorrelation function using statsmodels library"""

import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
import functions.initial as initial
import functions.metropolis as metropolis
import functions.cluster as cluster

width = 16  # width of the lattice
betaJ = 0.33

# Metropolis-Hastings ACF
print('MH mags')
mag_MH = []
MH_n_steps = 12800  # number of steps to run each simulation for Metropolis-Hastings algorithm
MH_burn_in = 12800  # euilibration time for Metropolis-Hastings algorithm

lattice = initial.create_lattice(width, 0)
metropolis.n_MH_moves(lattice, width, betaJ, MH_burn_in) # equilibration

for step in range(MH_n_steps):
    metropolis.MH_flip(lattice, width, betaJ)
    mag_MH.append(abs(initial.magnetisation(lattice)))
    if step%10 == 0:
        print(step)

#Wolf ACF
#print('Wolff mags')
#p_add = 1 - np.exp(-2*betaJ)
#mag_wolff = []
#wolff_n_steps = 500 # number of moves for wolff algorithm
#wolff_burn_in = 100 # equilibration time for wolff

#lattice = initial.create_lattice(width, 0)
#cluster.n_wolff_moves(lattice, p_add, wolff_burn_in) # equilibration

#for step in range(wolff_n_steps):
#    cluster.wolff_flip1(lattice, p_add)
#    mag_wolff.append(abs(initial.magnetisation(lattice)))
#    if step%10 == 0:
#        print(step)

ACF_MH = initial.autocorrelation(mag_MH)
ACF_MH2 = initial.autocorrelation_2(mag_MH,15360)
#ACF_wolff = initial.autocorrelation(mag_wolff)

#data = [ACF_MH,ACF_wolff]
#np.save('acf_samemoves4', data)


plt.figure()
plt.plot(ACF_MH2, label='autocorr2', color='blue')
plt.plot(ACF_MH, label='autocorr1(statsmodel)', color='red')
plt.axhline(y=np.exp(-1), color='black', linestyle='--', label='1/e')
plt.xlabel('Lag time')
plt.ylabel('ACF')
plt.title('ACF for MH vs Wolff at beta = {}'.format(betaJ))
plt.legend()
plt.show()