#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as pl
import sys

# Read in size of system
n = int(sys.argv[1])

# Initialise Hamiltonian matrix
H = np.zeros ( (n,n) )

# Read in coordinate file
coordfile = sys.argv[2]

print "Using file: " , coordfile

# Read in cell dimensions
cell=[float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5])]

print("Cell dimensions: ",cell)

# Read in LUMOs (assuming no energetic disorder)
LUMO=float(sys.argv[6])

# Load C60 locations from coordinate file. Format:-
#  X Y Z
#  Assuming angstroms.
locations=np.loadtxt(coordfile) 

locations=locations/cell # scale to fractional coordinates

distancematrix=locations[:,None,...]-locations[None,...] # rolled over
# Calculate distance matrix with Numpy functional programming methods.

PBCS=True
if (PBCS==True):
    distancematrix[distancematrix<0.5]+=1.0 #minimum image convention
    distancematrix[distancematrix>0.5]-=1.0 #minimum image convention

distancematrix*=cell # scale back to real coordinates
locations*=cell # scale from fractional coordinates to real distances

H=np.apply_along_axis(np.linalg.norm,2,distancematrix) # distances via linalg norm command on suitables axes
# elements in H are now euler distances between those sites {i,j}

# Get transfer integrals for off-diagonal elements
J0=10
BETA=0.6
H=J0*np.exp(-BETA*H) # calculate transfer integrals with isotropic exponential form

# Fill diagonal elements with site LUMOs
np.fill_diagonal(H,LUMO)

# Solve Hamiltonian for evals and evecs
evals,evecs=np.linalg.eigh(H)

# Find Inverse Participation Ratio for lowest state
p=evecs[:,0]*evecs[:,0]
PR=0
PR_all=0

for i in range(0,n):
	PR+=p[i]*p[i]

IPR=1./PR

print "Inverse participation ratio for lowest state is: ", IPR

IPR_all=np.zeros(n)

# Find IPR for all states and save
for i in range(0,n):
	p_all=evecs[:,i]*evecs[:,i]
	for j in range(0,n):
		PR_all+=p_all[j]*p_all[j]

	IPR_all[i]=1/PR_all
	PR_all=0


np.savetxt("IPR_%s.dat"%(coordfile),IPR_all,delimiter=' ',newline='\n')

# Save eigenvalues 
np.savetxt("Eigenvalues_%s.dat"%(coordfile),evals,delimiter=' ',newline='\n')








