#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
import random
from scipy.sparse.linalg import eigsh
import sys


def setup_ham(js,sigma,homo):
	"Setup Hamiltonian with site energies on the diagonal and coupling between nearest neighbours on the off-diagonals"

	ham = js
    
	np.fill_diagonal(ham,homo)
    
	if sigma!=0:
		np.fill_diagonal(ham,np.random.normal(loc=homo,scale=sigma,size=N))

	return ham


def plot_js(ham):
    "Plot J values on a 2D grid"

    np.fill_diagonal(ham, 0.0) 
    fig=pl.figure()
    pl.axes().set_aspect('equal') 
    pl.title("Off-diagonal elements of Hamiltonian")
    pl.imshow(ham,interpolation='nearest', cmap=pl.cm.PuBuGn
    pl.colorbar()
    pl.show()
    
    return ham


def SCpolarongenerator(ham,ham_p,site,alpha):
    "Generate a 'polaron' by allowing the site energies to be self-consistently perturbed in proportion to the charge density on the site, until convergence of lowest energy state"

    SCFSTEPS = 1000
    init_diagonal = np.diagonal(Ham)
    
    np.fill_diagonal(Ham_p,init_diagonal)

    Ham_p[site,site]-=alpha

    evals,evecs=eigsh(Ham_p,1,which='LM',tol=1E-4) 

    for i in range(SCFSTEPS): 
        polaron=evecs[:,0]*evecs[:,0] # Find lowest energy state charge density
        np.fill_diagonal(Ham_p,init_diagonal-alpha*polaron) # Deepen site energies in proportion to density
        pvals,pvecs=eigsh(Ham_p,1,which='LM',tol=1E-4)  # Resolve Hamiltonian
        if np.isclose(pvals[0],evals[0],rtol=1e-8,atol=1e-8):  # Repeat until convergence of lowest state energy
            break
        evals=pvals
        evecs=pvecs

    print "SCF steps to localise: ", i

    return ham_p


def solve_ham(ham):
    "Solve the Hamiltonian"
    evals,evecs = np.linalg.eigh(ham)

    return evals,evecs


def get_occs(evecs,state):
    "Find occupations of each site"
    p=evecs[:,state]*evecs[:,state]

    return p


def get_ipr(evecs,state,n):
    "Calculate the Inverse Participation Ratio, varying between 1 and N, depending on the number of sites the charge is delocalised over"
    p=evecs[:,state]*evecs[:,state]
    pr=0

    for i in range(0,n):
        pr+=p[i]*p[i]

    ipr=1./pr

    return ipr


def max_occupied_site(evecs,state):
    "Find the maximally occupied site"
    p=evecs[:,state]*evecs[:,state]
    max_occ_site=np.argmax(p)

    return max_occ_site


def dos(evals,nbins):
    "Get histogram of eigenvalues to form a dos"
    hist,bin_edges=np.histogram(evals,bins=nbins)

    fig = pl.figure()
    pl.plot(bin_edges[0:nbins],hist[0:nbins])

    pl.show()

    fig.savefig("Spiro_crystal_DOS.png")
    
    return bin_edges[0:nbins],hist[0:nbins]


def main():

	jif = np.loadtxt(sys.argv[1])

	n = int(sys.argv[2])

	print "N: ", n

	m = len(jif[:,2])

	js=np.zeros((n,n))             

	for i in range(0,m):
		js[int(jif[i,0]),int(jif[i,1])]=jif[i,2]
		js[int(jif[i,1]),int(jif[i,0])]=jif[i,2]


	alpha = 0.0		  # Polaron formation energy
	sigma = 0.0       # Energetic disorder
	homo = -5.0

	ham = setup_ham(js,sigma,homo)
	ham_p = ham+0.0

	print "Max J= ", np.max(ham)

	evals,evecs=solve_ham(ham)

	evals_all=[]
	iprs=[]

	evals,evecs = solve_ham(ham_p)
	print "IPR:", get_ipr(evecs,0,n)

	dos(evals,50)
	

main()



