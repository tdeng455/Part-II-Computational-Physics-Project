import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt

"""Creates a lattice and defines useful functions such as calculating neighbouring sites, 
calculating the sum of neighbouring spin sites and computes the total magnetisation of the lattice"""

# Creates the lattice
def create_lattice(width,type=0):
    lattice = np.array((width,width))   
    if type == 1:
        lattice = np.ones((width,width))
    elif type == 0:
        lattice = np.random.choice([1,-1],(width,width))
    elif type == -1:
        lattice = (-1)*np.ones((width,width))
    else:
        raise ValueError("Invalid type. Type should be 0, 1, or -1.")
    return lattice

# Plots the lattice with black as spin up, and white as spin down
def plot_lattice(lattice,ax,title):
    ax.matshow(lattice, vmin=-1, vmax=1, cmap=plt.cm.binary)
    ax.title.set_text(title)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_xticks([])

# Returns the coordinates of the neighbouring sites of a given coordinate with periodic boundrary conditions
def get_neighbouring_sites(i,j,width):
    neighbours = [((i-1)%width, j), ((i+1)%width, j), (i, (j-1)%width), (i, (j+1)%width)]

    return neighbours

# Returns the sum of the spins of the neighbouring sites of a given coordinate
def neighbouring_spins_sum(i,j,lattice,width):
    neighbours = get_neighbouring_sites(i,j,width)
    sum_neighbours = sum(lattice[i] for i in neighbours)
    return sum_neighbours

# Computes the overall magnetisation of the lattice
def magnetisation(lattice):
    Mag = lattice.sum()/(len(lattice)**2)
    return Mag