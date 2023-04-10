"""onsager.py
Plots the analytic solution for the mean magnetisation of a 2D Ising model with no external field
"""

import numpy as np
import matplotlib.pylab as plt

T = np.linspace(0.1,4,500)
T_crit = 2/np.log(1+np.sqrt(2))

def Mag_pos(x):
    if x >= T_crit:
        return 0
    else:
        return ((1 - np.sinh(2/x)**(-4))**(1/8))
def Mag_neg(x):
    if x >= T_crit:
        return 0
    else:
        return -1*((1 - np.sinh(2/x)**(-4))**(1/8))
  
M_pos = np.zeros(len(T))
M_neg = np.zeros(len(T))

for i in range(len(T)):
    M_pos[i] = Mag_pos(T[i])
    M_neg[i] = Mag_neg(T[i])

fig = plt.figure()
plt.title("Onsager solution for mean magentisation in the limit of large lattice")
plt.plot(T, M_pos, label=None, color='b')
plt.plot(T, M_neg, label=None, color='b')
plt.xlabel ('Temperature, T')
plt.ylabel ('Magnetisation, M')
plt.show()