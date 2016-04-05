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

# Open the log files
    molA_parser=ccopen("%s"%(MOLA_proj))
    molB_parser=ccopen("%s"%(MOLB_proj))
    molAB_parser=ccopen("%s"%(MOLAB_proj))

# Parse the relevant data
    molA=molA_parser.parse()
    molB=molB_parser.parse()
    molAB=molAB_parser.parse()

    print ("Parsed...")

# Size of basis sets
    nbasisA=molA.nbasis
    nbasisB=molB.nbasis
    
    print "nbasisA: ", nbasisA
    print "nbasisB: ", nbasisB
    
    Nbasis=nbasisA+nbasisB

# Position of HOMO    
    nhomoA=molA.homos
    nhomoB=molB.homos
    
    print "nhomoA: ", nhomoA
    print "nhomoB: ", nhomoB
    
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
    MOspairorth=np.dot(Dpair,MOsAB)

# Calculate the Fock matrix    
    B=np.dot(MOsorth.T,MOspairorth)
    Evals=np.diagflat(EvalsAB)    
    F=np.dot(np.dot(B,Evals),B.T)

# Output the HOMO-HOMO and LUMO-LUMO coupling elements fromt the Fock matrix
    print "HOMO-HOMO coupling: ", F[nhomoB+nbasisA,nhomoA]
    print "LUMO-LUMO coupling: ", F[nhomoB+nbasisA+1,nhomoA+1]


ProJ()
