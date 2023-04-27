"""Plots magnetisation against number of MH sweeps at different temps to 
   illustrate equilibration of the algorithm
"""

import numpy as np
rng = np.random.default_rng()
import functions.initial as initial
import functions.metropolis as metropolis

def equilibration_MH(width, betaJ, n_sweeps, filename):   

    lattice_types = [(0, "Random uniform"), (1, "All up"), (-1, "All down"), (2, "Anti-aligned")]
    magnetisation_data = []
    for i, (lattice_type, name) in enumerate(lattice_types):
        lattice = initial.create_lattice(width, type=lattice_type)
        print(i)
        
        #Measurement
        N = width**2
        mag = np.abs(initial.magnetisation(lattice))
        magnetisations = [mag]
        for j in range(n_sweeps*N):
            metropolis.MH_flip(lattice, width, betaJ)
            if j in np.arange(n_sweeps)*N:
                mag = np.abs(initial.magnetisation(lattice))
                magnetisations.append(mag)
        
        magnetisation_data.append([name, magnetisations])
        #ax.plot(range(n_sweeps+1), magnetisations, label=name)

    np.save(filename, np.array(magnetisation_data, dtype=object))
    #ax.set_xlabel("Number of Sweeps")
    #ax.set_ylabel("Absolute Magnetisation, |M|")
    #ax.set_title(r"Metropolis-Hastings Equilibration for $\beta$ = {}".format(betaJ))
    #ax.legend()

equilibration_MH(40,0.33,200, 'MH_equilibration_data_0.33')
equilibration_MH(40,1,200, 'MH_equilibration_data_1.00')