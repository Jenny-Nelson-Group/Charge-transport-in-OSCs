# Transform a molecule to its equivalent positions in a unit cell using equivalent positions from the space group, and translate to make a 3x3x3 supercell. 

Spacegroup.py defines a subset of the transformations from space groups (can easily be extended from http://homepage.univie.ac.at/nikos.pinotsis/spacegroup.html)

Transform_molecules.py takes as arguments:

1. An xyz file of one molecule
2. The cell lengths of the unit cell (which can be extracted form a cif file with Extract_cell_lengths.sh)
3. The space group (as defined in Spacegroup.py)
4. The file name for the final files

Outputs xyz files for all the transformed and translated molecules.




