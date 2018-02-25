#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
import random
from scipy.sparse.linalg import eigsh
import sys


def Setup_Ham(distance_matrix,sigma,lumo):
    "Setup Hamiltonian with site energies on the diagonal and coupling between nearest neighbours on the off-diagonals"
    
    H=np.apply_along_axis(np.linalg.norm,2,distancematrix)
    
    J0=10
    BETA=0.6
    H=J0*np.exp(-BETA*H)
    
    np.fill_diagonal(H,lumo)
    
    if sigma!=0:
        np.fill_diagonal(H,np.random.normal(loc=lumo,scale=sigma,size=N))

    return H


def SCpolarongenerator(Ham,Ham_p,site,alpha):
    "Generate a 'polaron' by allowing the site energies to be self-consistently perturbed in proportion to the charge density on the site, until convergence of lowest energy state"

    SCFSTEPS = 2000
    init_diagonal = np.diagonal(Ham)
    
    np.fill_diagonal(Ham_p,init_diagonal)

    Ham_p[site,site]-=alpha

    Evals,Evecs=eigsh(Ham_p,1,which='LM',tol=1E-4)  # Solve Hamiltonian


    for i in range(SCFSTEPS): # Number of SCF steps
        #print i
        polaron=Evecs[:,0]*Evecs[:,0] # Find lowest energy state charge density
        np.fill_diagonal(Ham_p,init_diagonal-alpha*polaron) # Deepen site energies in proportion to density
        pvals,pvecs=eigsh(Ham_p,1,which='LM',tol=1E-4)  # Resolve Hamiltonian
        if np.isclose(pvals[0],Evals[0],rtol=1e-8,atol=1e-8):  # Repeat until convergence of lowest state energy
            break
        Evals=pvals
        Evecs=pvecs

#print i

    return Ham_p



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



# Load data and apply PBCs

coordfile=sys.argv[1]
locations=np.loadtxt(coordfile)
N=len(locations)
cell=[float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4])]

locations=locations/cell

distancematrix=locations[:,None,...]-locations[None,...]

PBCS=False
if (PBCS==True):
    distancematrix[distancematrix<-0.5]+=1.0 #minimum image convention
    distancematrix[distancematrix>0.5]-=1.0 #minimum image convention

distancematrix*=cell # scale back to real coordinates
locations*=cell # scale from fractional coordinates to real distances


# Define parameters

ALPHA=0.006       # Energy of polaron formation
SIGMA=0.0       # Energetic disorder
LUMO=-3.7

# Main

H = Setup_Ham(distancematrix,SIGMA,LUMO)
Hp=H+0.0

EVALS=[]
IPRs=[]

for SITE in range(0,N):
    Ham_p = SCpolarongenerator(H,Hp,SITE,ALPHA)
    pvals,pvecs = Solve_Ham(Ham_p)
    EVALS=np.append(EVALS,pvals)
    ipr=IPR(pvecs,0)
    print SITE,ipr
    IPRs=np.append(IPRs,ipr)

OCCS=Get_occs(pvecs,0)

print IPRs

print "Average IPR: ", np.mean(IPRs)

np.savetxt("Evals_%s_%s_%s.dat"%(coordfile,ALPHA,SIGMA),EVALS)
np.savetxt("Evecs_%s_%s_%s.dat"%(coordfile,ALPHA,SIGMA),OCCS)
np.savetxt("IPRs_%s_%s_%s.dat"%(coordfile,ALPHA,SIGMA),IPRs)




