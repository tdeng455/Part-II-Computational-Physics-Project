"""Calculate and plot dynamic exponent via finite size scaling"""

import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
import functions.initial as initial
import functions.metropolis as metropolis
import functions.cluster as cluster
import functions.acfast as acfast

def compute_dynamic_exponent_MH(w_values, betaJ, burn_in, n_steps):
    
    #Generate data
    autotimes = []
    for L in w_values:
        mags=[]
        lattice = initial.create_lattice(L,1)
        metropolis.n_MH_moves(lattice, L, betaJ, burn_in)
        for i in range(n_steps):
            mags.append(initial.magnetisation(lattice))
            metropolis.MH_flip(lattice, L, betaJ)
            if i%100 == 0:
                mags.append(initial.magnetisation(lattice))
                print(i)
        autotime_L = acfast.autocorrelation_time(mags,len(mags))
        autotimes.append(autotime_L)
    print(autotimes)

    log_w = np.log(w_values)
    log_times = np.log(autotimes)
    coeffs = np.polyfit(log_w, log_times, deg=1)
    
    return log_w, log_times, coeffs

def compute_dynamic_exponent_wolff(w_values, betaJ, burn_in, n_steps):
    p_add = 1-np.exp(-2*betaJ)
    
    #Generate data
    autotimes = []
    for L in w_values:
        print(L)
        mags=[]
        lattice = initial.create_lattice(L,0)
        cluster.n_wolff_moves(lattice, p_add, burn_in)
        for i in range(n_steps):
            mags.append(initial.magnetisation(lattice))
            cluster.wolff_flip1(lattice, p_add)
            if i%1000 == 0:
                print(i)
        autotime_L = acfast.autocorrelation_time(mags,len(mags))
        autotimes.append(autotime_L)
    print(autotimes)

    log_w = np.log(w_values)
    log_times = np.log(autotimes)
    coeffs = np.polyfit(log_w, log_times, deg=1)

    return log_w, log_times, coeffs

log_w, log_times, coeffs = compute_dynamic_exponent_MH([10,20,30,40,50,60,70,100,150,200,500,1000],  0.440687, 50000, 100000)

data = [log_w,log_times]
np.save('log_width and log_times', data)

print(coeffs)
plt.figure()
plt.plot(log_w,log_times)
plt.show()