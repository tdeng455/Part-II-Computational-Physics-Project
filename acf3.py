"""Calculates and plots autocorrelation function using statsmodels library"""

import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
import functions.initial as initial
import functions.metropolis as metropolis
import functions.cluster as cluster

def acf_magnetisation_steps(width, betaJ, wolff=False, MH_n_steps = 500, 
                            MH_burn_sweeps = 200, wolff_n_steps = 200, 
                            wolff_burn_in = 1000):
    N = width**2

    if wolff == True:
        #Wolf ACF
        print('--------------Wolff case--------------')
        p_add = 1 - np.exp(-2*betaJ)
        mag_wolff = []

        print('Equilibrating...')
        lattice = initial.create_lattice(width, 0)
        cluster.n_wolff_moves(lattice, p_add, wolff_burn_in) # equilibration

        print('Generating Wolff series...')
        for step in range(wolff_n_steps):
            cluster.wolff_flip1(lattice, p_add)
            mag_wolff.append(abs(initial.magnetisation(lattice)))
            #if step%100 == 0:
            #    print(step)
        
        print('Calc of Wolff ACF...')
        ACF_wolff = initial.autocorrelation(mag_wolff)

        return ACF_wolff
    
    else:
        # Metropolis-Hastings ACF
        print('--------------Metropolis case--------------')
        mag_MH = []

        print('Equilibrating...')
        lattice = initial.create_lattice(width, 0)
        metropolis.n_MH_moves(lattice, width, betaJ, MH_burn_sweeps*N) # equilibration

        print('Generating MH series...')
        for step in range(MH_n_steps):
            metropolis.MH_flip(lattice, width, betaJ)
            mag_MH.append(abs(initial.magnetisation(lattice)))
            #if step%100 == 0:
            #    print(step)

        print('Calc of Metropolis-Hastings ACF...')
        ACF_MH = initial.autocorrelation(mag_MH)

        return ACF_MH

def acf_magnetisation_sweeps(width, betaJ, n_sweeps, MH_burn_sweeps = 500, 
                             wolff_burn_steps=2000, wolff=False):
    
    N = width**2

    if wolff == False: 
        
        #------------Metropolis ACF-------------#

        print('--------------Metropolis case--------------')
        mag_MH = []
        sweeps_data_MH = []
        MH_n_steps = n_sweeps*N  # number of steps to run each simulation for MH
        MH_burn_in = MH_burn_sweeps*N  # equilibration time for MH

        print('Equilibrating...')
        lattice = initial.create_lattice(width, 0)
        metropolis.n_MH_moves(lattice, width, betaJ, MH_burn_in) 

        print('Generating MH series...')
        for step in range(MH_n_steps):
    
            sweeps_data_MH.append(step/N)
            metropolis.MH_flip(lattice, width, betaJ)
            mag_MH.append(np.abs(initial.magnetisation(lattice)))

        print('Calculating MH ACF...')
        ACF_MH = initial.autocorrelation(mag_MH)
        #ACF_MH_2 = initial.autocorrelation_2(mag_MH, len(mag_MH))

        return ACF_MH, sweeps_data_MH
    
    else: 
        
        #------------Wolff ACF-------------#

        print('--------------Wolff case--------------')
        p_add = 1 - np.exp(-2*betaJ)

        #finding average cluster size
        print('Finding average cluster size...')
        lattice = initial.create_lattice(width, 0)
        cluster.n_wolff_moves(lattice, p_add, 5000)
        total_flips = cluster.n_wolff_moves(lattice, p_add, 5000)
        cluster_size = total_flips/5000

        mag_wolff = []
        sweeps_data_wolff = []
        wolff_n_steps = round(n_sweeps*N/cluster_size) # number of moves for wolff

        print('Equilibrating...')
        lattice = initial.create_lattice(width, 0)
        cluster.n_wolff_moves(lattice, p_add, wolff_burn_steps)

        print('Generating Wolff series...')
        num_flips = 0
        sweeps = 0
        for i in range(wolff_n_steps):
            num_flips += cluster.wolff_flip1(lattice, p_add)
            if num_flips > N:
                sweeps += num_flips/N
                #print(sweeps)
                sweeps_data_wolff.append(sweeps)
                m = np.abs(initial.magnetisation(lattice))
                mag_wolff.append(m)
                num_flips = 0

        print('Calculating Wolff ACF...')
        ACF_wolff = initial.autocorrelation(mag_wolff)
        #ACF_wolff_2 = initial.autocorrelation_2(mag_wolff, len(mag_wolff))

        return ACF_wolff, sweeps_data_wolff

"""
beta = 0.44
widths = [16,32,64,128,256,512]
ACF_MH_average = []
#metropolis
for w in widths:
    print('width = ', w)
    acf1 = acf_magnetisation_steps(w,beta,False)
    acf2 = acf_magnetisation_steps(w,beta,False)
    acf3 = acf_magnetisation_steps(w,beta,False)
    acf4 = acf_magnetisation_steps(w,beta,False)
    acf5 = acf_magnetisation_steps(w,beta,False)

    print('averaging')
    ACF_MH_average.append(np.mean(([acf1,acf2,acf3,acf4,acf5]), axis=0))

for i in range(len(ACF_MH_average)):
    plt.plot(ACF_MH_average[i], label=widths[i])
plt.legend()
plt.yscale("log")
plt.show()
"""