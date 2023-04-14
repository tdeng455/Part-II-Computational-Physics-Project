"""generates data for cluster size as a function of temperature"""

import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
import functions.initial as initial
import functions.cluster as cluster

width = 50
N = width**2
temps = np.linspace(1.5,3.5,13)
burn_in = 500
n_moves = 5000

cluster_sizes = []

for temp in temps:
    p_add = 1 - np.exp(-2/temp)
    lattice = initial.create_lattice(width, 0)
    cluster.n_wolff_moves(lattice, p_add, n_moves)
    total_flips = 0
    for i in range(n_moves):
        total_flips += cluster.n_wolff_moves(lattice,p_add,1)
cluster_sizes.append(total_flips/(n_moves))