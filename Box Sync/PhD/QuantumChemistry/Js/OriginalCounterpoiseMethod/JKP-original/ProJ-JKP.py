import numpy as np
import re

def TransFort(numFort):
	numSane=numFort[0:-4]+"E"+numFort[-3:len(numFort)]
	return float(numSane)

# NB: splitline modifies the list words that it is given
def splitline(line,n,words):
	if len(line) < n:
		return words
	else:
		numFort=line[0:n]
		
		words.append( float( TransFort( line[0:n] ) ) )
		splitline( line[n:len(line)], n, words)
		
def readorb(filename, n):
	psi=np.zeros([n,n])
	eps=np.zeros([n,n])

	f = open(filename, 'r')
	
	#skip the first line
	f.readline()
	
	#cycle for orbitals
	for i in range(0,n):
		orbnrg=f.readline()
	#	print "orbital nrg ",i, TransFort(orbnrg[-16:-1])
		eps[i,i] = TransFort(orbnrg[-16:-1])
		
		#cycle the basis set 
		j=0
		while j < n :
			line=f.readline()
			l=list()
			splitline(line,15, l)	
			
			for count in range (0, len(l)):
				psi[i, j] = l[count]
			#	print l[count], j
				j=j+1			

	f.close() 
	return [psi,eps]


def readS(filename,n):
	S=np.zeros([n,n])

	f = open(filename, 'r')

	line=()
	
	#skip all lines till the overlap matrix
	while line != " *** Overlap *** \n":
		line=f.readline()
	
	#outer loop over bs
	j=0
	imin=0	
	while j<n:

		#skip teh index at the beg
		line=f.readline()
		#inner loop
		i=imin
	#	print "start inner loop"
		while i<n:
			line=f.readline()
			purge=line[7:-1]
		#	print purge
			l=list()
			splitline(purge,14,l)
		#	print l
			
			#now add the elements to S
			for count in range(0,len(l)):
				S[i,j+count] = l[count]
				S[j+count,i] = l[count]	
			i=i+1
		j=j+5
		imin=imin+5
	f.close()
	return S 

#this is a little wasteful, could parse the file in one go. however it is simpler to just parse it several times
def GetNBasis(namefile):
	f = open(namefile, 'r')
	sNBasis = re.compile('NBasis')
	sNum = re.compile('[0-9]+')
	word=()
	while f:
		word=f.readline()
		test=sNBasis.search(word)
		if test != None:
			hits = sNum.findall(word)
			return int(hits[0])
			
	
	

def CalcJ(namejob):

	nbs=GetNBasis("dim/"+namejob+"dim.log")
	#call the aweful read functions
	[psipart1, nrgpart1]=readorb("part1/fort.7",nbs)
	[psipart2, nrgpart2]=readorb("part2/fort.7",nbs)
	[psidim, nrgdim]    =readorb("dim/fort.7",nbs)
	S = readS("dim/ethdim.log",nbs)

	#write the matrix of MOs including weighting from the S matrix:
	DimPro = np.transpose(np.dot(psidim, S))

	#project the orbitals of the two parts:
	Psi1DimBS = np.dot( psipart1, DimPro)
	Psi2DimBS = np.dot( psipart2, DimPro)

	#print nrgdim

	#determine the transfer integral
	JAB=np.dot(   np.dot(Psi1DimBS,nrgdim), np.transpose(Psi2DimBS) )
	JAA=np.dot(  np.dot(Psi1DimBS,nrgdim), np.transpose(Psi1DimBS) )
	JBB=np.dot(  np.dot(Psi2DimBS,nrgdim), np.transpose(Psi2DimBS) )
	return [JAB, JAA, JBB, Psi1DimBS, Psi2DimBS]

