#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
import scipy
from scipy import linalg, matrix
import sys



#Define null space function
def null(X,eps=1e-15):
    u, s, vh = scipy.linalg.svd(X)                      #Single value decomposition
    null_mask = (np.absolute(s)/np.max(s) <= eps)       #Check for values below tol
    null_space = scipy.compress(null_mask, vh, axis=0)  #Find corresponding vectors
    n,m=np.shape(null_space)
    for i in range(0,n):                                #Choose only solutions that
        positive=np.greater_equal(null_space[i,:],np.zeros(m))#can be probabilities
        negative=np.less_equal(null_space[i,:],np.zeros(m))
        if (np.all(positive)==True or np.all(negative)==True):
            Solution=np.absolute(null_space[i,:])
    Solution=Solution/np.sum(Solution)                  #Ensure sum(Pi)=1
    return scipy.transpose(Solution)



#Plot rate distribution to check
def PlotRates(X):
    logX=np.zeros((N,N))

    for i in range(0,N):
        for j in range(0,N):
            if X[i,j]!=0:
                logX[i,j]=np.log10(np.absolute(X[i,j]))
    logX=logX[np.nonzero(logX)]
    print logX
    hist_logX,bins=np.histogram(logX,bins=100)
    print bins
    fig=pl.figure()
    pl.bar(bins[0:99],hist_logX[0:99],width=0.01)

    pl.show()


#-------------------------------------------------------------------------------#


#Load data
N=int(sys.argv[1])
Jif=np.loadtxt(sys.argv[2])
coordfile=np.loadtxt(sys.argv[3])
#LUMOs=np.loadtxt(sys.argv[4])

#Define variables
A=np.zeros((N,N))
Lambda_inner=0.1
Lambda_outer=0.2
hbar=6.582*10**-16
e=1
kb=8.617*10**-5
T=300
F=[50000,0,0]       #In V/cm

Lambda=Lambda_inner+Lambda_outer

F_mag=np.linalg.norm(F)                 #In V/cm

print "F= ", F_mag, "V/cm"

#Enantiopure
#cell= [  [12.99,           0.0000000000,         0.0000000000,],
#       [0.0000000000,        17.88,              0.0000000000,],
#       [0.0000000000,        0.0000000000,         7.48] ]

#Racemic
cell= [  [11.2929,         0.0000000000,         0.0000000000,],
       [-0.524863,        10.691,         0.0000000000,],
       [2.23624,         2.1321,        6.97203] ]


#J_nonzero=Jif[:,2]
#print "Non-zero Js: ", J_nonzero[np.nonzero(J_nonzero)]

print "J file: ", Jif

Size=len(Jif[:,2])

print "Size of J file: ", Size

orderedJs=np.argsort(Jif[:,2])


print "Top Js: ", Jif[orderedJs[Size-5:Size],2], "at", Jif[orderedJs[Size-5:Size],0], Jif[orderedJs[Size-5:Size],1]

J=np.zeros((N+1,N+1))

for i in range(0,Size): 
    J[Jif[i,0],Jif[i,1]]=Jif[i,2]

#print "J: ", J


#PBCs

coordfile=coordfile*10**-8              #Convert from A to cm

coordfile=np.inner(np.linalg.inv(cell),coordfile).T   #scale to fractional coordinates
distancematrix=coordfile[:,None,...]-coordfile[None,...]

distancematrix[distancematrix<0.5]+=1.0 #minimum image convention
distancematrix[distancematrix>0.5]-=1.0

distancematrix=np.inner(distancematrix,cell).T # scale back to real coordinates
coordfile=np.inner(coordfile,cell).T

distancematrix=distancematrix.T

#print "Distance Matrix", distancematrix

dim=0

for i in range(0,3):
    if coordfile[1,i]!=0:
        dim+=1

#print "J", J[np.nonzero(J)] 

#PlotRates(J)

#Calculate matrix of rates
for i in range(0,N):
    for j in range(0,N):
        if (i!=j): #and np.absolute(Jif[i,j])>10e-5):
#            J[i,j]=J[j,i]
            deltaE=0
            d=np.absolute(distancematrix[i,j,:])
            field=np.dot(F,d)            
#print "exp= ", np.exp(-(((deltaE-field)+Lambda)**2)/(4*Lambda*kb*T))
            #print "J^2/hbar= ",(Jif[i,j]*Jif[i,j])/hbar
            A[i,j]=((J[i,j]*J[i,j])/hbar)*((np.pi/(Lambda*kb*T))**0.5)*np.exp(-(((deltaE-field)+Lambda)**2)/(4*Lambda*kb*T))
#	    if A[i,j]!=0:print A[i,j]

#print "Rate matrix: ", np.nonzero(A)

#print "Max Rate: ", np.max(A)

#PlotRates(A)

for i in range(0,N):
    A[i,i]=-np.sum(A[:,i])

print "Rates: ", A[np.nonzero(A)]

P=null(A)

print "P= ", P

print "AP= ", np.inner(A,P) 

r=np.zeros((N,N))

for i in range(0,N):
    for j in range(0,N):
        if i!=j:
            r[i,j]=np.linalg.norm(distancematrix[i,j,:])

#print "Distance: ", r

mu=0

for i in range(0,N):
    for j in range(0,N):
        if i!=j:
            mu+=np.dot(distancematrix[i,j,:],F)*P[i]*A[j,i]

mu=mu/F_mag


print "Mobility= ", mu, "cm^2 / Vs"

