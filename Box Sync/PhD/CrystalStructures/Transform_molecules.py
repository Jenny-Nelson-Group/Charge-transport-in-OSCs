#!/usr/bin/env python

# Transform molcules to equivalent positions and then translate, then make xyz files of each in supercell for calculation of J

import numpy as np
import sys
import Spacegroup


# --------------------------------- Read in data --------------------------------------- #


data = sys.argv[1]

# Extract raw data from xyz file (i.e. discard atom symbols)
LABELS = np.genfromtxt(data,usecols=0,dtype=str)
COORDFILE = np.genfromtxt(data)[:,1:]

# print "Coordinates: ", coordfile
N=int(sys.argv[2])

# Cell lengths
x=float(sys.argv[3])
y=float(sys.argv[4])
z=float(sys.argv[5])
 
# Read in space group: Pbca, Cc, P21c, P21, Pna21, Pbcn, P212121, C2c, P_1,Pca21, C2
SG = sys.argv[6]

# File name for output files
FILENAME = sys.argv[7]

# Cell matrix
CELL= [  [x,           0.0000000000,         0.0000000000,],
[0.0000000000,        y,              0.0000000000,],
[0.0000000000,        0.0000000000,         z] ]


# -------------------------------------------------------------------------------------- #

def transform_and_translate(coordfile,cell):

	# Transform to fractional coordinates
	coordfile=np.inner(np.linalg.inv(cell),coordfile).T

	#print "Fractional coordinates: ", coordfile

	#Add a column of ones for matrix multiplication as in Vesta
	coordfile=np.c_[coordfile,np.ones(N)]

	# Transform matrices

	Transforms = np.array(Spacegroup.choose(SG))

	print Transforms 

	# Translation matrices

	Translations=[ [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],
    	          [[1,0,0,1],[0,1,0,0],[0,0,1,0],[0,0,0,1]],
    	          [[1,0,0,0],[0,1,0,1],[0,0,1,0],[0,0,0,1]],
    	          [[1,0,0,0],[0,1,0,0],[0,0,1,1],[0,0,0,1]],
    	          [[1,0,0,-1],[0,1,0,0],[0,0,1,0],[0,0,0,1]],
    	          [[1,0,0,0],[0,1,0,-1],[0,0,1,0],[0,0,0,1]],
    	          [[1,0,0,0],[0,1,0,0],[0,0,1,-1],[0,0,0,1]],
    	          [[1,0,0,1],[0,1,0,1],[0,0,1,0],[0,0,0,1]],
    	          [[1,0,0,1],[0,1,0,-1],[0,0,1,0],[0,0,0,1]],
    	          [[1,0,0,-1],[0,1,0,1],[0,0,1,0],[0,0,0,1]],
    	          [[1,0,0,-1],[0,1,0,-1],[0,0,1,0],[0,0,0,1]],
    	          [[1,0,0,1],[0,1,0,0],[0,0,1,1],[0,0,0,1]],
    	          [[1,0,0,1],[0,1,0,0],[0,0,1,-1],[0,0,0,1]],
    	          [[1,0,0,-1],[0,1,0,0],[0,0,1,1],[0,0,0,1]],
    	          [[1,0,0,-1],[0,1,0,0],[0,0,1,-1],[0,0,0,1]],
    	          [[1,0,0,0],[0,1,0,1],[0,0,1,1],[0,0,0,1]],
    	          [[1,0,0,0],[0,1,0,1],[0,0,1,-1],[0,0,0,1]],
    	          [[1,0,0,0],[0,1,0,-1],[0,0,1,1],[0,0,0,1]],
    	          [[1,0,0,0],[0,1,0,-1],[0,0,1,-1],[0,0,0,1]],
    	          [[1,0,0,1],[0,1,0,-1],[0,0,1,-1],[0,0,0,1]],
    	          [[1,0,0,-1],[0,1,0,-1],[0,0,1,1],[0,0,0,1]],
    	          [[1,0,0,-1],[0,1,0,1],[0,0,1,-1],[0,0,0,1]],
    	          [[1,0,0,-1],[0,1,0,1],[0,0,1,1],[0,0,0,1]],
    	          [[1,0,0,1],[0,1,0,-1],[0,0,1,1],[0,0,0,1]],
    	          [[1,0,0,1],[0,1,0,1],[0,0,1,-1],[0,0,0,1]],
    	          [[1,0,0,1],[0,1,0,1],[0,0,1,1],[0,0,0,1]],
    	          [[1,0,0,-1],[0,1,0,-1],[0,0,1,-1],[0,0,0,1]]
    	          ]

	all_transformations = np.zeros((1,N,3))

	for displace_n,displace in enumerate(Translations):
	    for transform_n,transform in enumerate(Transforms):
	        # Transform coordinates
			trans=np.inner(transform,coordfile).T
	        # Translate in one of 6 directions
			trans=np.inner(displace,trans).T
	        # Remove column of ones
			trans=trans[:,[0,1,2]]
	        # Scale up to real coordinates
			trans=np.inner(trans,cell)
			all_transformations=np.vstack((all_transformations,trans[None,...]))

	print np.shape(all_transformations[1:])			

	return all_transformations[1:]


def check_in_cell(coordfile):
	coordfile_fractional=np.inner(np.linalg.inv(CELL),coordfile).T
	atoms=[]
	for i in range(0,N):
		if 0<coordfile_fractional[i,0]<1 and 0<coordfile_fractional[i,1]<1 and 0<coordfile_fractional[i,2]<1:
			atoms=np.append(atoms,i)

	return len(atoms)


def unit_cell(all_molecules):
	in_cell=[0,0]
	for i in range(0,len(all_molecules)):
		check=check_in_cell(all_molecules[i,:,:])
		if check>0:
			in_cell=np.vstack((in_cell,[i,check]))
	in_cell=in_cell[1:]

	sorted_idx=np.argsort(in_cell[:,1])
	in_cell_sorted=in_cell[sorted_idx][::-1]

	unitcell=[0,0]
	
	for i in range(len(in_cell_sorted)):
		for j in range(len(all_molecules)/27):
			if in_cell_sorted[i,1]%j not in unitcell[:,0]:
				unitcell=np.vstack((unitcell,[in_cell_sorted[i,1]%j,in_cell_sorted[i,1]]))

	print unitcell


	return unitcell


	




def save_files(file_to_save,labels,FILENAME,i):
	# Put labels back to save in xyz format
	trans=np.array(zip(labels,file_to_save[:,0],file_to_save[:,1],file_to_save[:,2]),dtype=		[('labels','S8'),('trans[:,0]',float),('trans[:,1]',float),('trans[:,2]',float)])
	# Save
	np.savetxt("%s_%d.xyz"%(FILENAME,i),trans,delimiter=" 	",fmt=["%s"]+["%f"]+["%f"]+["%f"])

all_mols = transform_and_translate(COORDFILE,CELL)
unit_cell(all_mols)


