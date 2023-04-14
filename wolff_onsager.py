"""Plots the analyitical and wolff simulated solutions for the average magnetisation of the lattice
"""

import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
import functions.initial as initial
import functions.cluster as cluster

betaJ_analytical = np.linspace(0,1.5,200)
beta_crit = 0.5*np.log(1+np.sqrt(2))

#defining analytical magnetisation solution
def Mag_pos(x):
    if x <= beta_crit:
        return 0
    else:
        return ((1 - np.sinh(2*x)**(-4))**(1/8))
def Mag_neg(x):
    if x <= beta_crit:
        return 0
    else:
        return -1*((1 - np.sinh(2*x)**(-4))**(1/8))

x = betaJ_analytical  
M_analytical_pos = np.zeros(len(x))
M_analytical_neg = np.zeros(len(x))

for i in range(len(x)):
    M_analytical_pos[i] = Mag_pos(x[i])
    M_analytical_neg[i] = Mag_neg(x[i])
    
avg_times = [75000,150000,225000,300000] # points to sample magnetisation

# different beta
width = 25
lattice = initial.create_lattice(width, type=0)
M_mcmc = []
betaJ_mcmc = np.linspace(0,1.5,60)
for i in betaJ_mcmc:
    print(i)
    p_add = 1-np.exp(-2*i)
    M_mcmc.append(cluster.compute_M_avg_wolff(lattice, p_add, avg_times=avg_times)) 

fig = plt.figure()
plt.plot(betaJ_analytical, M_analytical_pos, label='Analytical', color='b')
plt.plot(betaJ_analytical, M_analytical_neg, label=None, color='b')
plt.plot(betaJ_mcmc, M_mcmc, label='Wolff MC Method', color='black')
plt.xlabel (r'Inverse temperature $\beta$')
plt.ylabel ('Magnetisation, M')
plt.legend(loc='upper left')
plt.title(r'Simulated and anaytical solutions for $\langle M \rangle$ against $\beta$')
plt.savefig('MH_onsager.png')
plt.show()