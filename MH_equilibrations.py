"""Plots magnetisation against number of MH sweeps at different temps to 
   illustrate equilibration of the algorithm
"""

import numpy as np
rng = np.random.default_rng()
import matplotlib.pylab as plt
import functions.initial as initial
import functions.metropolis as metropolis

def plot_equilibration_MH(width, betaJ, n_sweeps, ax):   

    lattice_types = [(0, "Random uniform"), (1, "All up"), (-1, "All down"), (2, "Anti-aligned")]
    for i, (lattice_type, name) in enumerate(lattice_types):
        lattice = initial.create_lattice(width, type=lattice_type)

        # Measurement stage
        N = width**2
        mag = np.abs(initial.magnetisation(lattice))
        magnetisations = [mag]
        for j in range(n_sweeps*N):
            metropolis.MH_flip(lattice, width, betaJ)
            if j in np.arange(n_sweeps)*N:
                mag = np.abs(initial.magnetisation(lattice))
                magnetisations.append(mag)

        ax.plot(range(n_sweeps+1), magnetisations, label=name)

    ax.set_xlabel("Number of Sweeps")
    ax.set_ylabel("Absolute Magnetisation, |M|")
    ax.set_title(r"Metropolis-Hastings Equilibration for $\beta$ = {}".format(betaJ))
    ax.legend()

fig, axs = plt.subplots(1,2, figsize=(14,5))
plot_equilibration_MH(50,0.33,200,axs[0])
plot_equilibration_MH(50,1,200,axs[1])
plt.savefig('MH_equilibration_temps.png',dpi='figure')
plt.show()