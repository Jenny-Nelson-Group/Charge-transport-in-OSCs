#!/usr/bin/env python

# Calculate the mobility given a crystal structure and the electronic couplings (Js) from the central unit cell to all others in a 3x3x3 supercell. The J file should be in the format $i $j J, where $i and $j are the indices for each molecule, and this code is set up the take advantage of the symmetry in a crystal and duplicate the Js between additional unit cells. The codes Transform_molecules.py and Translate_molecules.py can be used to make up the unit cells, given symmetry operations for the crystal structure. 

# The code works by using the populated J matrix to find the corresponding rates in the prescence of an electric field (although field strength can be set to zero and the mobility calculated via the Einstein relation). 

# A master equation for rates is then set up and solved: AP=0, where A is the rate matrix containing the sum of all rates away from a molecule on the diagonal and the rate to each other molecule on the off diagonals. This is used to find the velocity, which can then be used to find the mobility via the field.

# The diffusion coefficient is calculated with a model based on a random walk, which is converted to a mobility with the Einstein relation.



import numpy as np
import matplotlib.pyplot as pl
import scipy
from scipy import linalg, matrix
from scipy.sparse.linalg import svds
import sys
from sympy import Matrix


# ----------------  Define null space function for solving AP=0 ------------------ # 

def null(X,eps=1e-12):
    Solution=0
    u, s, vh = scipy.linalg.svd(X)     # Single value decomposition
    null_mask = (np.absolute(s) <= eps)

    null_space = scipy.compress(null_mask, vh, axis=0)  # Find corresponding vectors

    n,m=np.shape(null_space)
    for i in range(0,n):           # Choose only solutions that can be probabilities
        positive=np.greater_equal(null_space[i,:],np.zeros(m))     
        negative=np.less_equal(null_space[i,:],np.zeros(m))
        if (np.all(positive)==True or np.all(negative)==True): 	   # Ensure sum(Pi)=1
            Solution=np.absolute(null_space[i,:])
            Solution=Solution/np.sum(Solution)
    return scipy.transpose(Solution)

# ------------------- Find matching rows (for populating full J matrix ------------ #


def find_rows(a, b):
    
    row_match = np.array(np.all((np.isclose(a[:,None,:],b[None,:,:],rtol=1e-3,atol=1e-10)),axis=-1).nonzero()).T.tolist()


    return row_match 


#-------------------------------- Main Code ------------------------------------- #


# Load data

N=int(sys.argv[1])                # Number of molecules
coordfile=np.loadtxt(sys.argv[2]) # read in coordinate file (in Angstroms)
Jif=np.loadtxt(sys.argv[3])		  # Read in J file (with columns $i $j J)
F_mag=float(sys.argv[4])          # Magnitude of field vector (in V/cm). Set to zero for no field.
theta=np.radians(float(sys.argv[5]))         # Angles to define direction of field
phi=np.radians(float(sys.argv[6]))

# Define constants

A=np.zeros((N,N))              # Initialise rate matrix
Lambda_inner=0.1               # Define inner and outer reorganisation energies
Lambda_outer=0.2
hbar=6.582*10**-16             # in eV.s
e=1                            # Charge on electron in eV/V 
kb=8.617*10**-5				   # in eV/K
T=300                          # in K
M=N/27      				   # M = number of molecules in unit cell (for 3x3x3 supercell)

#print "M= ", M

# Define cartesian coordinates of field vector (assuming from origin)

F=[F_mag*np.sin(phi)*np.cos(theta),F_mag*np.sin(phi)*np.sin(theta),F_mag*np.cos(phi)]                   

#print "Field: ", F

F_MAG=np.linalg.norm(F)        # Calculate magnitude of field (= F_mag)

Lambda=Lambda_inner+Lambda_outer  # Calculate total lambda


Size=len(Jif[:,2])

orderedJs=np.argsort(Jif[:,2])

#print "Top Js: ", Jif[orderedJs[Size-10:Size],2], "at", Jif[orderedJs[Size-10:Size],0], Jif[orderedJs[Size-10:Size],1]

J=np.zeros((N,N))              # Set up Js in matrix

for i in range(0,Size): 
    J[Jif[i,0],Jif[i,1]]=Jif[i,2]

# Convert from Angstroms to cm

#print "Non zero Js: ", np.count_nonzero(J)

coordfile=coordfile*10**-8

distancematrix=coordfile[:,None,...]-coordfile[None,...]

#print "Distance matrix: ", np.shape(distancematrix)

# Find sites in 1st unit cell (by finding rows of matrix containing more than a critical value of non-zeros)

unit_cell_sites=[]

for i in range(0,N):
    if len(J[np.nonzero(J[i,:])])>M+1:
    	unit_cell_sites=np.append(unit_cell_sites,[i])

unit_cell_sites = [int(x) for x in unit_cell_sites]      # Set elements as integers

unit_cell_vecs=np.zeros((N,N,3))
rest_distance_matrix=np.zeros((N,N,3))

unit_cell_vecs[0:M,:,:]=distancematrix[0:M,:,:]
rest_distance_matrix[M:N,:,:]=distancematrix[M:N,:,:]

for i in range(0,M):
	for j in range(M,N):
		J_match = find_rows(unit_cell_vecs[i,:,:],rest_distance_matrix[j,:,:])
		if unit_cell_vecs[i,J_match[0][0],:].all()!=0:
			for match in range(0,len(J_match)):
				if i!=match:
					J[j,J_match[match][1]]=J[i,J_match[match][0]]
					J[J_match[match][1],j]=J[J_match[match][0],i]
				#print "J= ", J[i,J_match[match,0]], "for ", unit_cell_vecs[i,J_match[match,0],:]

#print "Non-zero Js after filling: ", np.count_nonzero(J)


# Calculate matrix of rates

for i in range(0,N):
    for j in range(0,N):
        if (i!=j):
            deltaE=0
            d=distancematrix[i,j,:]
            field=np.dot(F,d)
            #print "field= ", field
            A[i,j]=((J[i,j]*J[i,j]))*((np.pi/(Lambda*kb*T))**0.5)*np.exp(-(((deltaE-field)+Lambda)**2)/(4*Lambda*kb*T))
            #print A[i,j]


for i in range(0,N):
	A[i,i]=-np.sum(A[:,i])

#print "A: ", A, "Non-zero: ", np.count_nonzero(A)


P_all=null(A)

#print "P: ", P_all


# Matrix of distances

r=np.apply_along_axis(np.linalg.norm,2,distancematrix)  


# -----Mobility from velocity (can only be calculated in prescence of a field ----- #

if F_MAG!=0:

	V=np.zeros(3)       # Find velcocity vector

	for i in range(0,N):
		for j in range(0,N):
			if i!=j:
				V+=distancematrix[i,j,:]*A[i,j]*P_all[i]


	V_F=np.dot(V,(F/F_MAG))       # V in direction of field

	mu = (V_F/F_MAG)/hbar                # Mobility in cm^2/Vs

	#print "At theta=",theta, ", phi=", phi , " Mobility= ", mu, "cm^2 / Vs"

	print theta, phi, mu


# ---------------------- Mobility from Einstein relation ------------------------- #

# Calculate D. Initialise.

D=0

# Unit vector to calculate mobility in the direction of

unit_vec = [np.sin(phi)*np.cos(theta),np.sin(phi)*np.sin(theta),np.cos(phi)]  

for i in range(0,N):
    for j in range(0,N):
        D+=0.5*P_all[i]*A[j,i]*(np.dot(distancematrix[j,i,:],unit_vec))**2

mu_ein= ((e*D)/(kb*T))/hbar
#print "At theta=",theta, ", phi=", phi, "Mobility from Einstein relation= ",mu_ein, "cm^2 / Vs"


