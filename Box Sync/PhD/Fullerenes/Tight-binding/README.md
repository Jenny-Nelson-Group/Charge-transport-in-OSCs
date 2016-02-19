# Set up Tight-binding Hamiltonian and solve
''For an input structure PDB file of coarse-grained fullerene assembly, this code exploits the isotropic coupling bewteen molecules to calculate the eigenvalues and eigenvectors though a tight-binding method: the code Evals.py outputs the Eigenvalues (use in conjunction with DOS.sh), while IPR.py calculates the size of the charge from the Eigenvectors with the Inverse Participation Ratio. The code Polaron_formation.py also includes polaron formation in a self-consistent manner. Written by Jarvist Frost and Beth Rice.''

The inputs for each code are the PDB structure file (where the number of molecules and cell dimensions are extracted by the ITIAM_final.sh script). The values for both energetic disorder (dx) and electron phonon coupling strength (alpha) can be defined in the script in eV.

Each code begins by reading in the structure data and converting it into a distance matrix, which is used to populate a tight-binding Hamiltonian, assuming that the electronic couplings scale exponentially with distance (as has been shown for fullerenes by JMF and others). In the simplest versions of the code (no polaron formation), this Hamiltonian is solved with the numpy linear algebra library. Evals.py outputs a list of Eigenvalues, which can be binned with Bin_data.py for plotting a histogram (the DOS). IPR.py outputs the inverse participation ratio for the lowest state and save the IPR for all states.

The Polaron_formation.py code uses a self-consistent method to 'trap' the charge, modelling the effect of a polaron. The code starts the charge density on one site and solves the TB Hamiltonian given this site deepened in energy by alpha. The code then recalculates the occupation of each site and deepens the energy of each site in proportion to the charge density on that site. This is repeated until convergence of the energy within a certain tolerance. 








