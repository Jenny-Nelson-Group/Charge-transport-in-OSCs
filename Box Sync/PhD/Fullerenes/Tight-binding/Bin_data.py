#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as pl
import sys


# Read in eigenvalues
evals = np.loadtxt(sys.argv[1])

# Read in adduct

adduct = sys.argv[2]

# Plot histogram

hist,bin_edges=np.histogram(evals,np.linspace(-0.11,0.36,100))

n=len(hist)

hist_data=np.array(zip(bin_edges[0:n],hist[0:n]),dtype=[('bin_edges[0:n]',float),('hist[0:n]',float)])

print len(evals)

np.savetxt("Binned_data_same_bins_%s.dat"%(adduct),hist_data,delimiter=' ',newline='\n')
