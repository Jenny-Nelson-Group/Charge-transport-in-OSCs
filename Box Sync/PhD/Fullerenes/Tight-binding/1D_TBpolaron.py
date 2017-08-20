#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
import random


def Setup_Ham(E,Edis,J0,Jdis):
    "Setup tridiagonal Hamiltonian with site energies on the diagonal and coupling between nearest neighbours on the off-diagonals"
    
    if Edis!=0:
        Es=np.random.normal(loc=E,scale=Edis,size=N)
    else:
        Es=np.ones(N)*E
    if Jdis!=0:
        Js=np.random.normal(loc=J0,scale=Jdis,size=N-1)
    else:
        Js=np.ones(N-1)*J0
    
    H = np.diag(Es,0) + np.diag(Js,-1) + np.diag(Js,1)

    return H


def SCpolarongenerator(Ham,Ham_p):
    "Generate a 'polaron' by allowing the site energies to be self-consistently perturbed in proportion to the charge density on the site, until convergence of lowest energy state"

    SCFSTEPS = 100
    init_diagonal = np.diagonal(Ham)
    
    print init_diagonal

    np.fill_diagonal(Ham_p,init_diagonal)

    Evals,Evecs=np.linalg.eigh(Ham_p)  # Solve Hamiltonian


    for i in range(SCFSTEPS): # Number of SCF steps
        print i
        polaron=Evecs[:,0]*Evecs[:,0] # Find lowest energy state charge density
        np.fill_diagonal(Ham_p,init_diagonal-ALPHA*polaron) # Deepen site energies in proportion to density
        pvals,pvecs=np.linalg.eigh(Ham_p) # Resolve Hamiltonian
        if np.isclose(pvals[0],Evals[0],rtol=1e-8,atol=1e-8):  # Repeat until convergence of lowest state energy
            break
        Evals=pvals
        Evecs=pvecs

    np.fill_diagonal(Ham_p,init_diagonal)

    return Ham_p,pvals,pvecs



def Solve_Ham(Ham):
    "Solve the Hamiltonian. Can be used even if no polaron formed"
    evals,evecs=np.linalg.eigh(Ham)

    return evals,evecs


def Get_occs(evecs,state):
    "Find occupations of each site"
    p=evecs[:,state]*evecs[:,state]

    return p

def IPR(evecs,state):
    "Calculate the Inverse Participation Ratio, varying between 1 and N, depending on the number of sites the charge is delocalised over"
    p=evecs[:,state]*evecs[:,state]
    PR=0

    for i in range(0,N):
        PR+=p[i]*p[i]

    IPR=1./PR

    return IPR


def Max_occupied_site(evecs,state):
    "Find the maximally occupied site"
    p=evecs[:,state]*evecs[:,state]
    Max=np.argmax(p)

    return Max


def DOS(evals,nbins):
    "Get histogram of eigenvalues to form a DOS"
    hist,bin_edges=np.histogram(evals,nbins=nbins)
    
    return bin_edges[0:nbins],hist[0:nbinss]


# Define parameters

ALPHA=0.2       # Energy of polaron formation
N=10            # Number of units in chain

# Main

H = Setup_Ham(-3,0.05,0.5,0.01)

Hp=H+0.0

Ham_p,pvals,pvecs = SCpolarongenerator(H,Hp)

evals,evecs=Solve_Ham(H)

energies = np.diagonal(H)

Pert_energies = np.diagonal(Ham_p)

Occs=Get_occs(evecs,0)

Pert_occs=Get_occs(pvecs,0)

print "IPR unperturbed: ", IPR(evecs,0)
print "IPR perturbed: ", IPR(pvecs,0)

fig=pl.figure()

pl.plot(np.linspace(0,N,N),Occs,"b")
pl.plot(np.linspace(0,N,N),energies,"b.")
pl.plot(np.linspace(0,N,N),Pert_occs,"r")
pl.plot(np.linspace(0,N,N),Pert_energies,"r.")

pl.show()
















