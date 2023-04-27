"""generates data for cluster size as a function of temperature"""

import numpy as np
rng = np.random.default_rng()
import functions.initial as initial
import functions.cluster as cluster

def av_cluster_sizes(lattice, p_add, burn_in, n_moves):
    
    cluster.n_wolff_moves(lattice, p_add, burn_in)
    
    total_flips = 0
    cluster_size = []
    for i in range(n_moves):
        total_flips += cluster.n_wolff_moves(lattice,p_add,1)
    cluster_size = (total_flips/(n_moves))

    return cluster_size

width = 40
N = width**2
temps = np.linspace(1.5,3.5,30)
burn_in = 200
n_moves = 500

cluster_sizes = []

for temp in temps:
    print(temp)
    p_add = 1 - np.exp(-2/temp)
    lattice = initial.create_lattice(width, 0)
    cluster_sizes.append(av_cluster_sizes(lattice, p_add, burn_in, n_moves))

data = [temps,cluster_sizes]
np.save('cluster_sizes_data', data)