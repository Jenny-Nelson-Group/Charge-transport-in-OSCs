#!/usr/bin/env python

import numpy as np
from scipy.sparse.linalg import eigsh
import sys


#--------------------------------------- Define functions ---------------------------------------------

# Fill the diagonal elements with site energy; for tight binding
def filldiagonal(Ham):
    np.fill_diagonal(Ham,-3.7)
    Ham_p=H+0.0 #no copy
    if dx!=0.0:np.fill_diagonal(Ham_p,np.random.normal(loc=-3.7,scale=dx,size=n))
    Ham=Ham_p+0.0
    Ham_diagonal=np.diagonal(Ham_p)
    np.savetxt("LUMOS_%s_%s_%s"%(coordfile,dx,ALPHA),Ham_diagonal,delimiter=' ',newline='\n')
    return Ham,Ham_p


#SC polaron self trapping effect with charge initially localised on each site

def SCpolarongenerator_allsites_sparse(Ham,Ham_p):
    SCFSTEPS = 1000
    S=SCend-SCstart 
    print "S= ", S
    siteEs=[]
    polarons=[]
    overlaps=[]
    pvecs_polaron=np.zeros((n,S))
    pvals_polaron=np.zeros((n,S))
    pvecs_size=np.zeros((n,n))
    evecs=np.zeros((n,n))
    Hamp_diag=np.zeros((n,S))
    init_diagonal = np.diagonal(Ham) #must use H0, otherwise is updated on every loop
    sites_occupied=[]   
    IPRs=[]

#    print init_diagonal

    for site in range(SCstart,SCend):
	np.fill_diagonal(Ham_p,init_diagonal) # reset Hamiltonian to unperturbed state
        Ham_p[site,site]-=ALPHA # apply local distortion
	Hp_diagonal=np.diagonal(Ham_p)
 
        # Adapted for Sparse eigensolver, only largest mag. eigenvalue
        # See: http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html#scipy.sparse.linalg.eigsh
        Evals,Evecs=eigsh(Ham_p,1,which='LM',tol=1E-4) # Default tol=0 = machine precision
        for i in range(SCFSTEPS): # Number of SCF steps
            polaron=Evecs[:,0]*Evecs[:,0] #lowest energy state electron density
#            print "sum(polaron): ",sum(polaron)
            np.fill_diagonal(Ham_p,init_diagonal-ALPHA*polaron)
            pvals,pvecs=eigsh(Ham_p,1,which='LM',tol=1E-4)
            if np.isclose(pvals[0],Evals[0],rtol=1e-8,atol=1e-8):
                break
            # Store these polaron values for the next loop comparison
            Evals=pvals
            Evecs=pvecs
   
        print i," SCFSTEPS to converge"
	print  "IPR of polaron ", state, "is ", IPR(pvecs,0)
        print "Site occupied is ", Max_occupied_site(pvecs,0)
        sites_occupied=np.append(sites_occupied,Max_occupied_site(pvecs,0))
        IPRs=np.append(IPRs,IPR(pvecs,0))

        PVALS,PVECS=np.linalg.eigh(Ham_p)    # Final calculation of all evals and evecs   
#        print "Evals: ", PVALS.T 
        pvals_polaron[:,site%S]=PVALS.T
#	print "pvals_polaron: ", pvals_polaron
        pvecs_polaron[:,site%S]=PVECS[:,0].T      #All lowest state polaron wavefuntions for calculation of J
        Hamp_diag[:,site%S]=np.diagonal(Ham_p)
        evecs=np.zeros((n,n))
    
	np.fill_diagonal(Ham_p,init_diagonal)

    

    np.savetxt("Evecs_%s_%s_%s_%s.dat"%(coordfile,ALPHA,dx,SCstart),pvecs_polaron,delimiter=' ',newline='\n')
    np.savetxt("Evals_%s_%s_%s_%s.dat"%(coordfile,ALPHA,dx,SCstart),pvals_polaron,delimiter=' ',newline='\n')
    np.savetxt("Perturbed_H_diagonal_%s_%s_%s_%s.dat"%(coordfile,ALPHA,dx,SCstart),Hamp_diag,delimiter=' ',newline='\n')
    np.savetxt("Polaron_Sites_%s_%s_%s_%s.dat"%(coordfile,ALPHA,dx,SCstart),sites_occupied,delimiter=' ',newline='\n')
    np.savetxt("IPR_%s_%s_%s_%s.dat"%(coordfile,ALPHA,dx,SCstart),IPRs,delimiter=' ',newline='\n')
    return Ham,Ham_p,pvals_polaron,pvecs_polaron


# Calculate IPR (inverse participation ratio)

def IPR(EVECS,STATE):
    p=EVECS[:,STATE]*EVECS[:,STATE]
    PR=0

    for i in range(0,n):
        PR+=p[i]*p[i]

    IPR=1./PR

    return IPR


# Solve final form of Hamiltonian
def solveHandHp(Ham):
    Evals,Evecs=np.linalg.eigh(Ham)
    
    return Evals,Evecs


#----------------------------------------- Read in inputs --------------------------------------------#

# Setup size of system to study...
# if present, read in number of sites from argument {1}
if len(sys.argv) > 1: n = int(sys.argv[1])
else: n=10

# Initialise our Hamiltonian matrix
H = np.zeros ( (n,n) )

# if present, read in coordinate filename from argument {2}
if len(sys.argv) > 2: coordfile = sys.argv[2]
else: coordfile="test.xyz"
print coordfile


# if present, read in cell coordinates from arguments {3,4,5}
if len(sys.argv) > 3: cell=[float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5])]
else: cell=[100,100,100]
#cell=[106.287,106.287,106.287] # Hard coded cell dimensions!!! FIXME
print("Cell dimensions: ",cell)

ALPHA = float(sys.argv[6])
dx = float(sys.argv[7])
state = float(sys.argv[8])
SCstart=int(sys.argv[9])
SCend=int(sys.argv[10])


# Load molecule locations from coordinate file. Format:-
#  X Y Z
#  Assuming angstroms.
locations=np.loadtxt(coordfile) #this defaults to reading them in as floats, which should be fine

locations=locations/cell # scale to fractional coordinates

distancematrix=locations[:,None,...]-locations[None,...] # rolled over
# Calculate distance matrix with Numpy functional programming methods.

# PBCs
PBCS=True
if (PBCS==True):
    distancematrix[distancematrix<0.5]+=1.0 #minimum image convention
    distancematrix[distancematrix>0.5]-=1.0 #minimum image convention

distancematrix*=cell # scale back to real coordinates
locations*=cell # scale from fractional coordinates to real distances

H=np.apply_along_axis(np.linalg.norm,2,distancematrix) # distances via linalg norm command on suitables axes
# elements in H are now euler distances between those sites {i,j}

J0=10
BETA=0.6
H=J0*np.exp(-BETA*H) # calculate transfer integrals with isotropic exponential form

Hp=H

# Initialise for outputs
evals=np.zeros(n)
pvals=np.zeros(n)
evecs=np.zeros((n,n))
pvecs=np.zeros((n,n))


#---------------------------Set up and solve, with polaron formation ----------------------------------------#

print "Generated Hamiltonian... "

H,Hp=filldiagonal(H)

print "Hamiltonian fully setup, time to solve!"


#Self-consistent polaron formation on every site
H,Hp,pvals_polaron,pvecs_polaron=SCpolarongenerator_allsites_sparse(H,Hp)

evals,evecs = solveHandHp(H)
print "Hamiltonian solved"


















