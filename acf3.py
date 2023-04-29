"""Calculates and plots autocorrelation function using statsmodels library"""

import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
import functions.initial as initial
import functions.metropolis as metropolis
import functions.cluster as cluster

def acf_magnetisation(width, betaJ, wolff=False, MH_n_steps = 200, MH_burn_in = 100000, 
                      wolff_n_steps = 200, wolff_burn_in = 100):
    
    if wolff == True:
        #Wolf ACF
        print('Wolff mags')
        p_add = 1 - np.exp(-2*betaJ)
        mag_wolff = []

        lattice = initial.create_lattice(width, 0)
        cluster.n_wolff_moves(lattice, p_add, wolff_burn_in) # equilibration

        for step in range(wolff_n_steps):
            cluster.wolff_flip1(lattice, p_add)
            mag_wolff.append(abs(initial.magnetisation(lattice)))
            if step%100 == 0:
                print(step)
        
        print('Calc of Wolff ACF')
        ACF_wolff = initial.autocorrelation(mag_wolff)

        return ACF_wolff
    
    else:
        # Metropolis-Hastings ACF
        print('MH mags')
        mag_MH = []

        lattice = initial.create_lattice(width, 0)
        metropolis.n_MH_moves(lattice, width, betaJ, MH_burn_in) # equilibration

        for step in range(MH_n_steps):
            metropolis.MH_flip(lattice, width, betaJ)
            mag_MH.append(abs(initial.magnetisation(lattice)))
            if step%100 == 0:
                print(step)

        print('Calc of Metropolis-Hastings ACF')
        ACF_MH = initial.autocorrelation(mag_MH)

        return ACF_MH


width = 16
betas = [0.01,0.33,0.44,0.55,1.00]
ACF_MH_average = []
#metropolis
for b in betas:
    print('beta = ', b)
    acf1 = acf_magnetisation(width,b,False)
    acf2 = acf_magnetisation(width,b,False)
    acf3 = acf_magnetisation(width,b,False)
    acf4 = acf_magnetisation(width,b,False)
    acf5 = acf_magnetisation(width,b,False)

    print('averaging')
    ACF_MH_average.append(np.mean(([acf1,acf2,acf3,acf4,acf5]), axis=0))

for i in range(len(ACF_MH_average)):
    plt.plot(ACF_MH_average[i], label=betas[i])
plt.legend()
plt.show()