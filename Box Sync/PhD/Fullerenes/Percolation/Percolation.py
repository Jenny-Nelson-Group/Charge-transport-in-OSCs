#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as pl
import sys

######################## Percolation/search functions ############################
############################  for finding clusters ###############################


# Define function to check if J>J_cutoff
def percolation(val,perc):
	if val > perc and val!=1:
		return 1
	else: 
		return 0

# Functions to search columns/rows of a matrix for J>J_cutoff, checking if alreading included in cluster (diagonal site element equal to 1)
def search_column(matrix,COLUMN):
	INDICES=[]
	for index,val in enumerate(matrix[:,COLUMN]):
		if percolation(val,J_cutoff)==1 and H[index,index]!=1:
			INDICES=np.append(INDICES,index)
	return INDICES

def search_row(matrix,ROW):
	INDICES=[]
	for index,val in enumerate(matrix[ROW,:]):
		if percolation(val,J_cutoff)==1 and H[index,index]!=1:
			INDICES=np.append(INDICES,index)
	return INDICES

# Function combining column and row searches to iteratively search column then row then column where J>J_cutoff. When true sets diagonal element for site to 1 to avoid going through same column/row multiple times
def search_all(ROW_INDICES,COLUMN_INDICES):
	print "Searching..."
	for ROW in np.unique(ROW_INDICES):
		COLUMN_INDICES=np.append(COLUMN_INDICES,search_row(H,ROW))
		H[ROW,ROW]=1
		for COLUMN in np.unique(COLUMN_INDICES):
			ROW_INDICES=np.append(ROW_INDICES,search_column(H,COLUMN))
			H[COLUMN,COLUMN]=1
	print "Done"
	return np.unique(ROW_INDICES), np.unique(COLUMN_INDICES),H


############################# Read in and load data ##############################

# Read in size of system
n = int(sys.argv[1])

# Initialise Hamiltonian matrix
H = np.zeros ( (n,n) )

# Read in coordinate file
coordfile = sys.argv[2]

print "Using file: " , coordfile

# Read in cell dimensions
cell=[float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5])]

# Define J cutoff

J_cutoff=float(sys.argv[6])

print("Cell dimensions: ",cell)

# Load C60 locations from coordinate file. Format:-
#  X Y Z
#  Assuming angstroms.
locations=np.loadtxt(coordfile) 


############################## Set up Matrix #####################################


locations=locations/cell # scale to fractional coordinates

distancematrix=locations[:,None,...]-locations[None,...] # rolled over
# Calculate distance matrix with Numpy functional programming methods.

PBCS=True
if (PBCS==True):
    distancematrix[distancematrix<0.5]+=1.0 #minimum image convention
    distancematrix[distancematrix>0.5]-=1.0 #minimum image convention

distancematrix*=cell # scale back to real coordinates
locations*=cell # scale from fractional coordinates to real distances

H=np.apply_along_axis(np.linalg.norm,2,distancematrix) # distances via linalg norm command on suitables axes
# elements in H are now euler distances between those sites {i,j}

# Get transfer integrals for off-diagonal elements
J0=10
BETA=0.6
H=J0*np.exp(-BETA*H) # calculate transfer integrals with isotropic exponential form

np.fill_diagonal(H,0)


################################ Find clusters ####################################

column_indices=[]
row_indices=[]
row_indices_old=[0]

# Search through a chosen column to find rows to search. Continue until all sites in cluster found
for i in range(0,1):
    row_indices=search_column(H,i)
    while len(np.unique(row_indices))!=len(np.unique(row_indices_old)):  
		row_indices_old=row_indices
		#print np.unique(row_indices_old)
		row_indices,column_indices,H=search_all(np.unique(row_indices),np.unique(column_indices))
		#print len(np.unique(row_indices))
		#print len(np.unique(column_indices))
    print "for column ", i, " cluster is: ", np.unique(row_indices), "size is: ", len(np.unique(row_indices))

				
			
			

















