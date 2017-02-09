#!/usr/bin/env python

# Transform molcules to equivalent positions and then translate, then make xyz files of each in supercell for calculation of J

import numpy as np
import sys
import Spacegroup


def read_in_files():
	data = sys.argv[1]

	# Extract raw data from xyz file (i.e. discard atom symbols)
	labels = np.genfromtxt(data,usecols=0,dtype=str)
	coordfile = np.genfromtxt(data)[:,1:]

	# print "Coordinates: ", coordfile
	n=int(sys.argv[2])

	# Cell lengths
	x=float(sys.argv[3])
	y=float(sys.argv[4])
	z=float(sys.argv[5])

	# Read in space group: Pbca, Cc, P21c, P21, Pna21, Pbcn, P212121, C2c, P_1, 		Pca21, C2
	SG = sys.argv[6]

	# File name for output files
	filename = sys.argv[7]



def transform_and_translate():

	# Cell matrix
	cell= [  [x,           0.0000000000,         0.0000000000,],
    	   [0.0000000000,        y,              0.0000000000,],
    	   [0.0000000000,        0.0000000000,         z] ]


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

	return trans








def save_files(file_to_save,labels,filename):
		
	# Put labels back to save in xyz format
	trans=np.array(zip(labels,file_to_save[:,0],file_to_save[:,1],file_to_save[:,2]),dtype=		[('labels','S8'),('trans[:,0]',float),('trans[:,1]',float),('trans[:,2]',float)])
	# Save
	np.savetxt("%s_%d.xyz"%(filename,transform_n+displace_n*m),trans,delimiter=" 	",fmt=["%s"]+["%f"]+["%f"]+["%f"])




