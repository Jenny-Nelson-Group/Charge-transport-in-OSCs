#!/usr/bin/env python

# Calculate the mobility given a crystal structure and the electronic couplings (Js) from the central unit cell to all others in a 3x3x3 supercell. The J file should be in the format $i $j J, where $i and $j are the indices for each molecule, and this code is set up the take advantage of the symmetry in a crystal and duplicate the Js between additional unit cells. The code Transform_molecules.py and can be used to make up the unit cells, given symmetry operations for the crystal structure.

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


# ----------------------  Define null space function for solving AP=0 --------------------------------- #

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


# -------------------------- Find matching rows (for populating full J matrix ------------------------- #

def find_rows(a,b):
    
    row_match = np.array(np.all((np.isclose(a[:,None,:],b[None,:,:],rtol=1e-3,atol=1e-10)),axis=-1).nonzero()).T.tolist()


    return row_match 


#----------------------------------  Find sites in 1st unit cell  ------------------------------------- #

def fill_Js(J):
    
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

    return J


# --------------------------------- Calculate matrix of rates ------------------------------------------ #

def Marcus_rates(J,field):
    
    A=np.zeros((N,N))

    for i in range(0,N):
        for j in range(0,N):
            if (i!=j):
                deltaE=0
                d=distancematrix[i,j,:]
                Field_d=np.dot(field,d)
                #print "field= ", field
                A[i,j]=((J[i,j]*J[i,j]))*((np.pi/(Lambda*kb*T))**0.5)*np.exp(-(((deltaE-Field_d)+Lambda)**2)/(4*Lambda*kb*T))
    return A



# --------------------------------- Solve Master equation to find Ps ----------------------------------- #


def Master_eq(A):

    for i in range(0,N):
        A[i,i]=-np.sum(A[:,i])

    P_all=null(A)

    #print "P: ", P_all

    return P_all



# ----------------- Mobility from velocity (can only be calculated in prescence of a field) ------------ #


def Mob(A,P_all,field_unit):

	V=np.zeros(3)       # Find velocity vector

	for i in range(0,N):
		for j in range(0,N):
			if i!=j:
				V+=distancematrix[i,j,:]*A[i,j]*P_all[i]


	V_F=np.dot(V,field_unit)       # V in direction of field

	mu = (V_F/F_mag)/hbar                # Mobility in cm^2/Vs
    
	print "Mobility= ", mu, " cm^2/Vs"

	return mu


# ------------------------------------ Mobility from Einstein relation --------------------------------- #

def Mob_einstein(A,P_all,theta,plane):

    # Calculate D. Initialise.
    D=0

    # Unit vector to calculate mobility in the direction of

    if plane=='ab':
        unit_vec = [np.cos(theta),np.sin(theta),0]
    elif plane=='bc':
        unit_vec = [0,np.cos(theta),np.sin(theta)]
    elif plane=='ac':
        unit_vec = [np.cos(theta),0,np.sin(theta)]


    for i in range(0,N):
        for j in range(0,N):
            D+=0.5*P_all[i]*A[j,i]*(np.dot(distancematrix[j,i,:],unit_vec))**2

    mu_ein= ((e*D)/(kb*T))/hbar

    print "Mobility from Einstein relation= ", mu_ein, " cm^2/Vs"


    return mu_ein

# -------------------------------------- Plot mobility vs angle  --------------------------------------- #


def Plot_mobility(Angles,Mobilities,plane):

    print "average: ", np.mean(np.absolute(Mobilities))

    print "min: ", np.min(np.absolute(Mobilities))

    print "max: ", np.max(np.absolute(Mobilities))

    fig=pl.figure()

    ax=pl.subplot(111, polar=True)
    ax.plot(Angles,Mobilities,'o',label='Mobility (cm$^2$V$^{-1}$s$^{-1}$)')
    ax.grid(True)
    #ax.set_rmax(0.00525)
    #ax.set_ylim(-2700,2700)
    #ax.set_yticks(np.arange(-2700,2700,1000))
    ax.legend(loc='upper left', bbox_to_anchor=(-0.2,1.1))

    pl.show()
    fig.savefig("%s_mobility_%s.pdf"%(filename,plane))


#------------------------------  Read in data/ set constants and variables ----------------------------- #


# Load data

coordfile=np.loadtxt(sys.argv[1]) # read in coordinate file (in Angstroms)
Jif=np.loadtxt(sys.argv[2])		  # Read in J file (with columns $i $j J)
F_mag=float(sys.argv[3])          # Magnitude of field vector (in V/cm). Set to zero for no field.
filename=sys.argv[4]
N=len(coordfile)                  # Number of molecules

# Define constants

A=np.zeros((N,N))              # Initialise rate matrix
Lambda_inner=0.1               # Define inner and outer reorganisation energies
Lambda_outer=0.2
hbar=6.582*10**-16             # in eV.s
e=1                            # Charge on electron in eV/V
kb=8.617*10**-5				   # in eV/K
T=300                          # in K
M=N/27      				   # M = number of molecules in unit cell (for 3x3x3 supercell)

print M

Lambda=Lambda_inner+Lambda_outer  # Calculate total lambda

Size=len(Jif[:,2])

orderedJs=np.argsort(Jif[:,2])

#print "Top Js: ", Jif[orderedJs[Size-10:Size],2], "at", Jif[orderedJs[Size-10:Size],0], Jif[orderedJs[Size-10:Size],1]

Js=np.zeros((N,N))              # Set up Js in matrix

for i in range(0,Size):
    Js[int(Jif[i,0]),int(Jif[i,1])]=Jif[i,2]

coordfile=coordfile*10**-8     # Convert from Angstroms to cm

distancematrix=coordfile[:,None,...]-coordfile[None,...]

J_all=fill_Js(Js)

Mob_ab=[]
Mob_bc=[]
Mob_ac=[]
angles=[]


for theta in np.linspace(0,2*np.pi,36):
    FIELDab=[F_mag*np.cos(theta),F_mag*np.sin(theta),0]
    FIELDbc=[0,F_mag*np.cos(theta),F_mag*np.sin(theta)]
    FIELDac=[F_mag*np.cos(theta),0,F_mag*np.sin(theta)]
    
    FIELDab_unit=[np.cos(theta),np.sin(theta),0]
    FIELDbc_unit=[0,np.cos(theta),np.sin(theta)]
    FIELDac_unit=[np.cos(theta),0,np.sin(theta)]
    
    Ratesab=Marcus_rates(J_all,FIELDab)
    Ratesbc=Marcus_rates(J_all,FIELDbc)
    Ratesac=Marcus_rates(J_all,FIELDac)
    
    Pab=Master_eq(Ratesab)
    Pbc=Master_eq(Ratesbc)
    Pac=Master_eq(Ratesac)
    
    mobilityab=Mob(Ratesab,Pab,FIELDab_unit)
    mobilitybc=Mob(Ratesbc,Pbc,FIELDbc_unit)
    mobilityac=Mob(Ratesac,Pac,FIELDac_unit)

    Mob_ab=np.append(Mob_ab,mobilityab)
    Mob_bc=np.append(Mob_bc,mobilitybc)
    Mob_ac=np.append(Mob_ac,mobilityac)

    angles=np.append(angles,theta)


Plot_mobility(angles,Mob_ab,'ab')
Plot_mobility(angles,Mob_bc,'bc')
Plot_mobility(angles,Mob_ac,'ac')













