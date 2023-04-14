"""Defines ACF and AC Time using faster methods than in 'initial.py' using statsmodels library"""

import numpy as np
from statsmodels.tsa.stattools import acovf #,acf

def autocorrelation(data, t_max):
    if t_max<len(data):
        N=t_max
    else:
        N=None
    autocov = acovf(data,adjusted=False,demean=True,fft=False,nlag=N)
    return autocov/autocov[0]

def autocorrelation_time(data, t_max):
    autocorr = autocorrelation(data, t_max)
    crit = np.exp(-1)
    time = np.argmin(autocorr>crit, axis=0)
    return time