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

# Read in space group: Pbca, Cc, P21c, P21, Pna21, Pbcn, P212121, C2c, P_1, 		Pca21, C2
SG = sys.argv[6]

# File name for output files
FILENAME = sys.argv[7]

# Cell matrix
CELL= [  [x,           0.0000000000,         0.0000000000,],
[0.0000000000,        y,              0.0000000000,],
[0.0000000000,        0.0000000000,         z] ]


# -------------------------------------------------------------------------------------- #

def transform_and_translate(coordfile):

	# Transform to fractional coordinates
	coordfile=np.inner(np.linalg.inv(cell),coordfile).T

	#print "Fractional coordinates: ", coordfile

	#Add a column of ones for matrix multiplication as in Vesta
	coordfile=np.c_[coordfile,np.ones(n)]

	# Transform matrices

	Transforms = Spacegroup.choose(SG)

	m=len(Transforms)

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

	all_transformations = ((1,N,3))

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
			all_transformations=np.vstack(all_transformations,trans[None,...])

	print all_transformations			

	return all_transformations


def check_in_cell(coordfile):
	coordfile_fractional=np.inner(np.linalg.inv(cell),coordfile).T

	for i in range(0,N):
    if 0<coordfile[i,0]<1 and 0<coordfile[i,1]<1 and 0<coordfile[i,2]<1:
        atoms=np.append(atoms,i)

	return len(atoms)


def unit_cell():
	





def save_files(file_to_save,labels,FILENAME,i):
	# Put labels back to save in xyz format
	trans=np.array(zip(labels,file_to_save[:,0],file_to_save[:,1],file_to_save[:,2]),dtype=		[('labels','S8'),('trans[:,0]',float),('trans[:,1]',float),('trans[:,2]',float)])
	# Save
	np.savetxt("%s_%d.xyz"%(FILENAME,i),trans,delimiter=" 	",fmt=["%s"]+["%f"]+["%f"]+["%f"])




