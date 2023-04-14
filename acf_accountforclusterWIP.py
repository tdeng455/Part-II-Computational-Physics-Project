import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
import functions.initial as initial
import functions.metropolis as metropolis
import functions.cluster as cluster
import functions.acfast as acfast

width = 40  # width of the lattice
n_sweeps = 100

# Metropolis-Hastings ACF
mag_MH = []
MH_n_steps = n_sweeps*width**2  # number of steps to run each simulation for Metropolis-Hastings algorithm
MH_burn_in = 20000  # equilibration time for Metropolis-Hastings algorithm
betaJ = np.log(1+np.sqrt(2))/2

lattice = initial.create_lattice(width, 0)
metropolis.n_MH_moves(lattice, width, betaJ, MH_burn_in) 
for step in range(MH_n_steps):

    metropolis.n_MH_moves(lattice, width, betaJ, 1)

    if step%width**2 == 0:
        print(step)
        m = np.abs(initial.magnetisation(lattice))
        mag_MH.append(m)

ACF_MH = acfast.autocorrelation(mag_MH, len(mag_MH))

#Wolf ACF
p_add = 1 - np.exp(-2*betaJ)
lattice = initial.create_lattice(width, 0)
cluster.n_wolff_moves(lattice, p_add, 500)
total_flips = 0
for i in range(5000):
    total_flips += cluster.n_wolff_moves(lattice,p_add,1)
cluster_size = (total_flips/(5000))
div = round(width**2/cluster_size)

mag_wolff = []
wolff_n_steps = round(n_sweeps*width**2/cluster_size) # number of moves for wolff algorithm
wolff_burn_in = 2000 # equilibration time for wolff

lattice = initial.create_lattice(width, 0)
cluster.n_wolff_moves(lattice, p_add, wolff_burn_in) # equilibration

for step in range(wolff_n_steps):
    cluster.n_wolff_moves(lattice, p_add, 1)  
    if step%div == 0:
        print(step)
        m = np.abs(initial.magnetisation(lattice))
        mag_wolff.append(m)

ACF_wolff = acfast.autocorrelation(mag_wolff, len(mag_wolff))

data=[ACF_MH,ACF_wolff]
np.save('acf_accountforsweeps', data)

plt.figure()
plt.plot(ACF_wolff, label='Wolff', color='blue')
plt.plot(ACF_MH, label='Metropolis-Hastings', color='red')
plt.xlabel('Sweeps of the latttice')
plt.ylabel('ACF')
plt.title(r'ACF for MH vs Wolff at critical temperature')
plt.legend()
plt.show()