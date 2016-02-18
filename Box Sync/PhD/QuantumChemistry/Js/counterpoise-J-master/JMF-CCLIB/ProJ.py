import numpy as np
from cclib.parser import ccopen

def CalcJ(namejob):
    molA_parser=ccopen("part1/"+namejob+"part1.log")
    molB_parser=ccopen("part2/"+namejob+"part2.log")
    molAB_parser=ccopen("dim/"+namejob+"dim.log")

    molA=molA_parser.parse()
    molB=molB_parser.parse()
    molAB=molAB_parser.parse()
    print ("Parsed...")

    nbs=molAB.nbasis
    if molA.nbasis!=molB.nbasis:
	print("Count of basis functions doesn't match. Failing.")
	return False

    for mole in molA,molB,molAB:
        if len(mole.atomcoords)!=1:
	    print(mole+" calculation appears to be an optimisation! Failing.")

#Take mocoeffs[0] - the alpha electrons to get matrix in correct order
    MolAB_Pro = np.transpose(np.dot(molAB.mocoeffs[0],molAB.aooverlaps))
#    print "DimPro via JKP:", DimPro, "DimPro2 via JMF", DimPro2

#    print "Psi1DimBS via JKP: ", Psi1DimBS
    PsiA_DimBS = np.dot(molA.mocoeffs[0], MolAB_Pro)
    PsiB_DimBS = np.dot(molB.mocoeffs[0], MolAB_Pro)
#    print "Dim Mocoeffs: ", molAB.moenergies[0]/27.211


#Note: moenergies in eV, so converted to Hartree for checking with JKP code
    JAB=np.dot(np.dot(np.diagflat(molAB.moenergies[0]/27.211),PsiA_DimBS), np.transpose(PsiB_DimBS) )
    JAA=np.dot(np.dot(np.diagflat(molAB.moenergies[0]/27.211),PsiA_DimBS), np.transpose(PsiA_DimBS) )
    JBB=np.dot(np.dot(np.diagflat(molAB.moenergies[0]/27.211),PsiB_DimBS), np.transpose(PsiB_DimBS) )

    print "JAB", JAB
    print "JAA", JAA
    print "JBB", JBB

    return [JAB, JAA, JBB]
