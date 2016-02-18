# Calculate percolation threshold of a system
'' Input a given system where electonic couplings (Js) are known between molecules. Calculates the percolation cluster size for different values of J_cutoff and therefore a percolation threshold in terms of J'

The Python code Percolation.py works for an input assembly with isotropic couplings (i.e. coupling only dependent on distance). The particular form of this function in the code is relevant for an assembly of course-grained fullerene molecules. An example pdb file is given. 

The code sets up a matrix of couplings between molecules based on this relationship and searches through each row and column, When J>J_cutoff, this molecule is added to the cluster and this is repeated until all clusters are found.

The script Percolation.sh tests different values of J_cutoff and puts the results in an output file for reading, where the relevant value is the size of the largest cluster. The script averages over many different assemblies. The script assumes you are in a directory that contains directories with the structure files in along with the code.

Finally, the script Extract_perc.sh extracts the value of J_cutoff along with the percolation cluster size for plotting.

