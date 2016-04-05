# Calculate the transfer integral J between two molcules with the counterpoise method 

Required

- Python
- cclib library for Python (install from here: http://cclib.sourceforge.net/wiki/index.php/Install)
- Gaussian
- xyz files of molA and molB


1. Run Gaussian calculations for molA, molB and molAB with HEADER_Js at the top of the com file. For molA and molB use ghost atoms (with -Bq after the atom symbol) of the other molecule to ensure the basis set is the same for all runs. This is for the b3lyp level of theory with basis set 6-31g*. IOp(6/7=3) ensures all molecular orbitals are printed. 

2. Copy the resultant log files into the same directory as Projective_method.py

3. Run 

python Counterpoise_method.py MolA.log MolB.log MolAB.log 


----------------------------------------------------------------------------------------------------------

This is based on work by James Kirkpatrick and Jarvist Frost.

The counterpoise method involves using ghost atoms to ensure molA, molB and molAB calculations are all done in the same basis. The non-orthogonalised orbitals of molA and molB are then projected onto the orbitals of molAB, giving phi_A and phi_B. The J matrix is then computed with <phi_A|E|phi_B>, where E is a diagonal matrix containing the eigenvalues of molAB





