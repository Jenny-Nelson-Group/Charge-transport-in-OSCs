import numpy as np
import scipy as sp
from scipy import linalg
import sys
from cclib.parser import ccopen


def ProJ():
# Read in molecule log files for projective method. Requires iop(3/33=1,6/7=3) in Gaussian header for calculation on each molecule + the pair
    MOLA_proj=sys.argv[1]
    MOLB_proj=sys.argv[2]
    MOLAB_proj=sys.argv[3]
    Degeneracy_HOMO=int(sys.argv[4])
    Degeneracy_LUMO=int(sys.argv[5])    # =0 for non-degenerate, 2,3 etc for doubly, triply etc

# Open the log files
    molA_parser=ccopen("%s"%(MOLA_proj))
    molB_parser=ccopen("%s"%(MOLB_proj))
    molAB_parser=ccopen("%s"%(MOLAB_proj))

# Parse the relevant data
    molA=molA_parser.parse()
    molB=molB_parser.parse()
    molAB=molAB_parser.parse()

    #print ("Parsed...")

# Size of basis sets
    nbasisA=molA.nbasis
    nbasisB=molB.nbasis
    
    #print "nbasisA: ", nbasisA
    #print "nbasisB: ", nbasisB
    
    Nbasis=nbasisA+nbasisB

# Position of HOMO    
    nhomoA=molA.homos
    nhomoB=molB.homos
    
    #print "nhomoA: ", nhomoA
    #print "nhomoB: ", nhomoB
    
# Get molecular orbitals. Need the transpose for our purposes.
    MOsA=(molA.mocoeffs[0]).T
    MOsB=(molB.mocoeffs[0]).T
    MOsAB=(molAB.mocoeffs[0]).T

# Get eigenvalues of pair    
    EvalsAB=molAB.moenergies[0]
    
# Get overlaps. These are symmetric so transpose not required in this case
    SA=molA.aooverlaps
    SB=molB.aooverlaps
    SAB=molAB.aooverlaps
 
# Set up matrices for MOs and S   
    MOs=np.zeros((Nbasis,Nbasis))
    S=np.zeros((Nbasis,Nbasis))
    
    MOs[0:nbasisA,0:nbasisA]=MOsA
    MOs[nbasisA:Nbasis,nbasisA:Nbasis]=MOsB

    S[0:nbasisA,0:nbasisA]=SA
    S[nbasisA:Nbasis,nbasisA:Nbasis]=SB

# Calculate upper diagonal matrix D, such that S=D.T*D for Lowdin orthogonalisation.       
    D=sp.linalg.cholesky(S)

    Dpair=sp.linalg.cholesky(SAB)

# Orthogonalise MOs matrix and MOsAB matrix
    MOsorth=np.dot(D,MOs)

    #print np.shape(MOsorth)

    MOspairorth=np.dot(Dpair,MOsAB)

    #print np.shape(MOspairorth)

# Calculate the Fock matrix    
    B=np.dot(MOsorth.T,MOspairorth)
    Evals=np.diagflat(EvalsAB)    
    F=np.dot(np.dot(B,Evals),B.T)


# Output the HOMO-HOMO and LUMO-LUMO coupling elements fromt the Fock matrix

    if Degeneracy_HOMO==0:
        print "HOMO-HOMO coupling: ", F[nhomoB+nbasisA,nhomoA]

    if Degeneracy_LUMO==0:
        print "LUMO-LUMO coupling: ", F[nhomoB+nbasisA+1,nhomoA+1]

# Degeneracies

    if Degeneracy_HOMO!=0:
        F_deg_HOMO=F[nhomoB+nbasisA-Degeneracy_HOMO+1:nhomoB+nbasisA+1,nhomoA-Degeneracy_HOMO:nhomoA]
        print "HOMO-HOMO coupling", (np.sum(np.absolute(F_deg_HOMO**2)))/Degeneracy_HOMO**2

    if Degeneracy_LUMO!=0:
        F_deg_LUMO=F[nhomoB+nbasisA:nhomoB+nbasisA+Degeneracy_LUMO,nhomoA:nhomoA+Degeneracy_LUMO]
        print "LUMO-LUMO coupling", (np.sum(np.absolute(F_deg_LUMO**2)))/Degeneracy_LUMO**2
		




ProJ()
