# Calculate transfer integral J between two molecules with the projective method

Required:

- Python 
- cclib library for Python (install from here: http://cclib.sourceforge.net/wiki/index.php/Install)
- Gaussian
- xyz files of molA and molB

1. Run Gaussian calculations for molA, molB and molAB with HEADER_Js at the top of the com file. This is for the b3lyp level of theory with basis set 6-31g*. IOp(3/33=1,6/7=3) ensures one-electron integrals and all molecular orbitals are printed. 

2. Copy the resultant log files into the same directory as Projective_method.py

3. Run 

python Projective_method.py MolA.log MolB.log MolAB.log 


-------------------------------------------------------------------------------------------------------
This is based on work by James Kirkpatrick and assistance from Jarvist Frost.

The projective method involves projecting the orbitals of the pair of molecules into the basis set of the individual molecules (for the case of orthogonal orbitals, as with the ZINDO method, more detail is given here http://arxiv.org/pdf/physics/0610288v1.pdf)

When the orbitals are not orthogonal, Lowdin orthogonalisation can be used to orthogonalise them.

 


 


