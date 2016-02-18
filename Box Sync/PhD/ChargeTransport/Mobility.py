#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
import scipy
from scipy import linalg, matrix
import sys



#Define null space function
def null(X,eps=1e-5):
    Solution=0
    u, s, vh = scipy.linalg.svd(X)                      #Single value decomposition
    null_mask = (np.absolute(s) <= eps)
    #print null_mask
    if not all(null_mask)==False:
        #print "Solution not found"
        return np.zeros(M)
    null_space = scipy.compress(null_mask, vh, axis=0)  #Find corresponding vectors
    n,m=np.shape(null_space)
#print n,m
    for i in range(0,n):                                #Choose only solutions that
        positive=np.greater_equal(null_space[i,:],np.zeros(m))#can be probabilities
        negative=np.less_equal(null_space[i,:],np.zeros(m))
        if (np.all(positive)==True or np.all(negative)==True): 	   #Ensure sum(Pi)=1
            #print null_space[i,:]
            Solution=np.absolute(null_space[i,:])
            Solution=Solution/np.sum(Solution)
    return scipy.transpose(Solution)


#Plot rate distribution to check
def PlotRates(X):
	logX=np.zeros((N,N))
	for i in range(0,N):
		for j in range(0,N):
			if X[i,j]!=0:
				logX[i,j]=np.log10(np.absolute(X[i,j]))
	logX=logX[np.nonzero(logX)]
    #print logX
	hist_logX,bins=np.histogram(logX,bins=100)
    #print bins
	fig=pl.figure()
	pl.bar(bins[0:99],hist_logX[0:99],width=0.01)
	pl.show()

    #fig.savefig("%s_dist.pdf"%(J_name))


#-------------------------------------------------------------------------------#


#Load data
N=int(sys.argv[1])
J_name=sys.argv[2]
Jif=np.loadtxt(sys.argv[2])
coordfile=np.loadtxt(sys.argv[3])
a=float(sys.argv[4])
b=float(sys.argv[5])
c=float(sys.argv[6])
F_mag=float(sys.argv[7])
f1=float(sys.argv[8])*F_mag
f2=float(sys.argv[9])*F_mag
f3=float(sys.argv[10])*F_mag
angle=float(sys.argv[11])

#Define variables
A=np.zeros((N,N))
Lambda_inner=0.1
Lambda_outer=0.2
hbar=6.582*10**-16
e=1
kb=8.617*10**-5
T=300
M=N/27      # m is number of molecules in unit cell

#print "m= ", M

F=[f1,f2,f3]       #In V/cm

#F=[0,0,0]

F_MAG=np.linalg.norm(F)

#print "F_MAG is: ", F_MAG

Lambda=Lambda_inner+Lambda_outer


#print "F= ", F, "V/cm"

cell= [ [a, 0.0, 0.0],
        [0.0, b, 0.0],
        [0.0, 0.0, c] ]


#Enantiopure exp
#cell= [  [16.542,         0.0000000000,         0.0000000000,],
 #      [0.0000000000,        15.035,         0.0000000000,],
 #      [-5.15855,         0.0000000000,        13.1326] ]



#print "J file: ", Jif

Size=len(Jif[:,2])

#print "Size of J file: ", Size

orderedJs=np.argsort(Jif[:,2])

#print "Top Js: ", Jif[orderedJs[Size-10:Size],2], "at", Jif[orderedJs[Size-10:Size],0], Jif[orderedJs[Size-10:Size],1]

J=np.zeros((N,N))

for i in range(0,Size): 
    J[Jif[i,0],Jif[i,1]]=Jif[i,2]

#print np.unique(J)
#print len(np.unique(J))


#PBCs

PBCs=False

#print coordfile

if PBCs==True:
	coordfile=np.inner(np.linalg.inv(cell),coordfile).T   #scale to fractional coordinates
	#print coordfile
	distancematrix=coordfile[:,None,...]-coordfile[None,...]
	distancematrix[distancematrix<0.5]+=1.0 #minimum image convention
	distancematrix[distancematrix>0.5]-=1.0
	distancematrix=np.inner(distancematrix,cell).T # scale back to real coordinates
	coordfile=np.inner(coordfile,cell).T
	distancematrix=distancematrix.T
else: distancematrix=coordfile[:,None,...]-coordfile[None,...] 

#print "Distance Matrix", distancematrix

#print "J", J

#PlotRates(J)

#Calculate matrix of rates
for i in range(0,N):
    for j in range(0,N):
        if (i!=j): #and np.absolute(Jif[i,j])>10e-5):
            #J[i,j]=J[j,i]
            deltaE=0
            d=distancematrix[i,j,:]
            #print "d= ", d
            field=np.dot(F,d)
            #print "field= ", field
            A[i,j]=((J[i,j]*J[i,j]))*((np.pi/(Lambda*kb*T))**0.5)*np.exp(-(((deltaE-field)+Lambda)**2)/(4*Lambda*kb*T))
#print A[i,j]

#print "J: ", len(J[np.nonzero(J)])

#print "A: ", len(A[np.nonzero(A)])

# Find sites in 1st unit cell (by finding rows of matrix containing more than a critical value of non-zeros)

unit_cell_sites=[]

for i in range(0,N):
    if len(A[np.nonzero(A[i,:])])>M:
        unit_cell_sites=np.append(unit_cell_sites,[i])

unit_cell_sites = [int(x) for x in unit_cell_sites]           # Set elements as integers

#unit_cell_sites=[0,1,2,3]

#print unit_cell_sites

# Set up matrix to solve for Ps

Mat=np.zeros((M,M))

for i in unit_cell_sites:
    Mat[i%M,i%M]=-np.sum(A[i,:])     # Fill diagonal with sum of rates away from molecule
    for j in range(0,N):
        Mat[i%M,j%M]+=A[j,i]

#P=np.zeros(N)

#P.fill(1./N)

#print A_cell

P= null(Mat)
#print "P= ", P

P_all=np.zeros(N)

for i in range(0,N):
    P_all[i]=P[i%M]/(N/4)


#print "P_all: ", P_all

#tile_P=np.tile(P,(N+1)/4) 
#P=tile_P/((N+1)/4)
#print P
#print sum(P_all)

#print "AP= ", np.inner(A,P)

r=np.zeros((N,N))

for i in range(0,N):
    for j in range(0,N):
        if i!=j:
            r[i,j]=np.linalg.norm(distancematrix[i,j,:])

#print "Distance: ", r

V=np.zeros(3)

for i in range(0,N):
    for j in range(0,N):
        if i!=j:
	    V+=distancematrix[i,j,:]*A[i,j]*P_all[j]

#print "v= ", V

mu=0

#for i in range(0,N):
#    for j in range(0,N):
#        if i!=j:
#            mu+=P_all[i]*A[j][i]*np.dot(distancematrix[j,i,:],(F/F_MAG))

#mu=mu/F_MAG

#print "mu= ", mu*10**-16

#Calculate D
#D=0
#Unit vector to calculate mobility on the direction of
#vec=F
#unit_vec=vec/np.linalg.norm(vec)

#for i in range(0,N):
#    for j in range(0,i):
#        if i!=j:
#            D+=0.5*P[i]*A[j,i]*(np.dot(distancematrix[j,i,:],unit_vec))**2


V_F=np.dot(V,(F/F_MAG))       #V in direction of field

print angle,((V_F/F_MAG)*10**-16)/hbar

#print "angle: ", angle , "mobility: ", (V_F/F_MAG)*10**-16

#print "Mobility= ", mu, "cm^2 / Vs"

#print "Mobility from Einstein relation= ", (e*D)/(kb*T) , "cm^2 / Vs"


