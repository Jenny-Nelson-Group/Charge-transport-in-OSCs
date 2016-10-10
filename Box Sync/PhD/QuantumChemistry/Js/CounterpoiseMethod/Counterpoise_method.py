import numpy as np
import scipy as sp
from scipy import linalg
import sys
from cclib.parser import ccopen


def CountPJ():
# Read in molecule log files for counterpoise method. Requires IOp(6/7=3) in Gaussian header + ghost atoms to make up basis sets in individual com files  
	MOLA_CP=sys.argv[1]
	MOLB_CP=sys.argv[2]
	MOLAB_CP=sys.argv[3]

# Open the log files
	molA_parser=ccopen("%s"%(MOLA_CP))
	molB_parser=ccopen("%s"%(MOLB_CP))
	molAB_parser=ccopen("%s"%(MOLAB_CP))

# Parse the relevant data
	molA=molA_parser.parse()
	molB=molB_parser.parse()
	molAB=molAB_parser.parse()

	print ("Parsed...")

# Size of basis sets
	Nbasis=molAB.nbasis

# HOMO index
	nhomoA=molA.homos
	nhomoB=molB.homos
	nhomoAB=molAB.homos

	nlumoA=nhomoA+1
	nlumoB=nhomoB+1

	print "HOMO A: ", nhomoA
	print "HOMO B: ", nhomoB
	print "HOMO AB: ", nhomoAB

# Every basis set should have the same size (the size of the pair) 
	if molA.nbasis!=molB.nbasis:
		print("Count of basis functions doesn't match. Failing.")
		return False

# Get molecular orbitals
	MOsA=molA.mocoeffs[0]
	MOsB=molB.mocoeffs[0]
	MOsAB=molAB.mocoeffs[0]

# Get eigenvalues of pair
	EvalsAB=molAB.moenergies[0]

	print "Energies: ", EvalsAB
	print "Gap: ", EvalsAB[nhomoAB+1]-EvalsAB[nhomoAB]
    
# Find HOMO and LUMO from energy splitting in dimer

	print "ESID HOMO-HOMO coupling", 0.5*(EvalsAB[nhomoAB]-EvalsAB[nhomoAB-1])
	print "ESID LUMO-LUMO coupling", 0.5*(EvalsAB[nhomoAB+2]-EvalsAB[nhomoAB+1])

# Calculate the molecular orbitals of A and B in the AB basis set

	SAB=molAB.aooverlaps

	MolAB_Pro = (np.dot(MOsAB,SAB)).T
	PsiA_AB_BS = np.dot(MOsA, MolAB_Pro)
	PsiB_AB_BS = np.dot(MOsB, MolAB_Pro)

# Calculate the matrix of transfer integrals

	JAB=np.dot(np.dot(PsiB_AB_BS,np.diagflat(EvalsAB)),PsiA_AB_BS.T)
	JAA=np.dot(np.dot(PsiA_AB_BS,np.diagflat(EvalsAB)),PsiA_AB_BS.T)
	JBB=np.dot(np.dot(PsiB_AB_BS,np.diagflat(EvalsAB)),PsiB_AB_BS.T)

	S = np.dot(PsiB_AB_BS,PsiA_AB_BS.T)

# Symmetric Lowdin transformation for required J

	J_eff_HOMO = (JAB[nhomoA,nhomoB]- 0.5*(JAA[nhomoA,nhomoA]+JBB[nhomoB,nhomoB])*S[nhomoA,nhomoB])/(1.0-S[nhomoA,nhomoB]*S[nhomoA,nhomoB])
	J_eff_LUMO = (JAB[nlumoA,nlumoB]- 0.5*(JAA[nlumoA,nlumoA]+JBB[nlumoB,nlumoB])*S[nlumoA,nlumoB])/(1.0-S[nlumoA,nlumoB]*S[nlumoA,nlumoB])


# Print the HOMO-HOMO and LUMO-LUMO coupling
	print "HOMO-HOMO coupling: ", J_eff_HOMO
	print "LUMO-LUMO coupling: ", J_eff_LUMO

	SA=molA.aooverlaps
	SB=molB.aooverlaps

	DA=sp.linalg.cholesky(SA)
	DB=sp.linalg.cholesky(SB)
	DAB=sp.linalg.cholesky(SAB)

	MOsA_orth = np.dot(DA,MOsA.T)
	MOsB_orth = np.dot(DB,MOsB.T)
	MOsAB_orth = np.dot(DAB,MOsAB.T)

	MolAB_Pro_orth=(np.dot(MOsAB_orth,SAB)).T
	PsiA_AB_BS_orth=np.dot(MOsA_orth,MolAB_Pro_orth)
	PsiB_AB_BS_orth=np.dot(MOsB_orth,MolAB_Pro_orth)

	JAB_lowd=np.dot(np.dot(PsiB_AB_BS_orth,np.diagflat(EvalsAB)),PsiA_AB_BS_orth.T)

	print "HOMO-HOMO coupling with Lowdin", JAB_lowd[nhomoA,nhomoB]
	print "LUMO-LUMO coupling with Lowdin", JAB_lowd[nlumoA,nlumoB]

CountPJ()



