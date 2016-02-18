#!/usr/bin/env python

# DBTW.py - generate DoS from CG Js read in from file
# Down by the water // My lovely daughter // I took her home

# Import our numeric library
import numpy as np
# Matplotlib
import matplotlib.pyplot as pl
# Sys for arg passing
import sys

import datetime # current date for log files etc.

# from IPython import embed# we do this so we can drop to interactive python for debugging; major Python coolio
 #  # --> embed() <-- just add this to any point in code, and TADA!

### Matplotlib setup
#Pretty colours; via http://blog.olgabotvinnik.com/post/58941062205/prettyplotlib-painlessly-create-beautiful-matplotlib
try:
    import brewer2mpl #'pip install brewer2mpl'
# Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
    colours = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
except ImportError: #If no brewer2mpl library
    #Otherwise, boring built in ones...
    colours='brgcmkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk' # Aaah, we fade to grey (fade to grey)
    print "Hey - no brewer2mpl (try 'pip install brewer2mpl'). Thus a simple palette."

# Matplotlib - initialise figure
fig=pl.figure()
pl.axes().set_aspect('equal') # Square data .'. square figure please

# Setup size of system to study...
if len(sys.argv) > 1: n = int(sys.argv[1])
else: n=10

# Initialise our Hamiltonian matrix
H = np.zeros ( (n,n) )

if len(sys.argv) > 2: edgesfile = sys.argv[2]
else: edgesfile="test.edges"

positionsfile = sys.argv[3]

# Load off-diagonal elements from text file. Format:-
#  (index_i,int) (index_j,int) (value,double) (value,double)
filein_edge=np.loadtxt(edgesfile,
                  dtype=[('f1', '<u4'), ('f2', '<u4'), ('f3',np.float64),('f4',np.float64)])  #specify datatypes: 4byte int, 4byte int, float64, float64

filein_pos=np.loadtxt(positionsfile,
                      dtype=[('f5', '<u4'),('f6',np.float64),('f7',np.float64),('f8',np.float64)])


for datum in filein_edge:
    #print datum
    H[datum[0],datum[1]]=datum[2] # Populate Hamiltonian with off diagonal elements
    H[datum[1],datum[0]]=datum[2]  #  Hermition...


print "Loaded Hamiltonian... "

#pl.title("Off-diagonal elements of Hamiltonian")
#pl.imshow(H,interpolation='nearest', cmap=pl.cm.PuBuGn) # 2D colourmap of Hamiltonian, nearest interpolation.
#pl.colorbar()
#pl.show()

# Fill the diagonal elements with site energy; for tight binding
#np.fill_diagonal(H, -6.0)


for j in range(0,n):
    H[j,j]=filein_edge[j][3]
#print H


print "Hamiltonian fully setup, time to solve!"
# OK; here we go - let's solve that TB Hamiltonian!
evals,evecs=np.linalg.eigh(H)


#Calculating centre of electrostatic density

centre_coords = np.zeros(3)
prob = np.zeros(n)
r=np.zeros(n)
distance_charge=np.zeros((n,2))
cum_prob=np.zeros(n)
sorted_charge=np.zeros(n)
polaron_size=np.zeros(n)
num=0.0

for i in range(0,100):

    for j in range(0,n):
        prob[j]=evecs[j,i]*evecs[j,i]
        centre_coords[0]+=filein_pos[j][1]*prob[j]
        centre_coords[1]+=filein_pos[j][2]*prob[j]
        centre_coords[2]+=filein_pos[j][3]*prob[j]


#Calculating how many molecules polarons localised over
        if prob[j]>0.01:num+=1
    #    print num



    for j in range(0,n):
        r[j]=np.sqrt((centre_coords[0]-filein_pos[j][1])*(centre_coords[0]-filein_pos[j][1])+(centre_coords[1]-filein_pos[j][2])*(centre_coords[1]-filein_pos[j][2])+(centre_coords[2]-filein_pos[j][3])*(centre_coords[2]-filein_pos[j][3]))
#    print r

    idx = np.argsort(r)
    sorted_r = r[idx]
    sorted_charge = prob[idx]

    cum_prob = np.cumsum(sorted_charge)
#    print cum_prob

    for j in range (0,n):
       if cum_prob[j]>0.95:break

    polaron_size[i] = sorted_r[j]

    centre_coords = np.zeros(3)     # Reset centre_coords
    num=0.0                         #Rest num

#print "Effective size of polaron for eigenvector 0 is", polaron_size[0]
#print "Effective size of polaron for eigenvector 5 is", polaron_size[5]
#print "Effective size of polaron for eigenvector 10 is", polaron_size[10]
#print "Effective size of polaron for eigenvector 50 is", polaron_size[50]





#Calculating electron probabilities
#max_coords = np.zeros(3)
#r=np.zeros(1000)

#for k in range(0,10):
    #Finding position of site with maxiumum probability
#    for i in range(0,n):
#        prob[i]=evecs[i,k]*evecs[i,k]
#    max_index=np.argmax(prob)
#    max_coords[0]=filein_pos[max_index][1]
#    max_coords[1]=filein_pos[max_index][2]
#    max_coords[2]=filein_pos[max_index][3]

#print "Prob=", prob[0]
#print "Max index=", max_index
#print "Max coordinates", max_coords


#Calculating an effective radius of the polaron by adding up weighted radii from centre found above
#    for j in range(0,n):
#        r[k]+=np.sqrt((max_coords[0]-filein_pos[j][1])*(max_coords[0]-filein_pos[j][1])+(max_coords[1]-filein_pos[j][2])*(max_coords[1]-filein_pos[j][2])+(max_coords[2]-filein_pos[j][3])*(max_coords[2]-filein_pos[j][3]))*prob[j]

#print r

#write eigenvalue vs size of polaron to files

#polaron_evals=np.zeros((n,2))
#for i in range (0,n):
#   polaron_evals[i,0]=evals[i]
#   polaron_evals[i,1]=r[i]

#np.savetxt("Polaron_vs_evals.dat",polaron_evals,delimiter=' ',newline='\n')

#Plot polaron size vs eigenvalue

fig1, ax=pl.subplots()
rects1 = pl.bar(evals[0:20],polaron_size[0:20],0.00001, color='b')
#rects2 = pl.bar(evals[0:20]+0.0001,r[0:20],0.0001, color='b')
pl.title("Size of polaron vs eigenvalue")
pl.xlabel("Eigenvalues (eV)")
pl.ylabel("Effective size of polaron ")
pl.ylim(0,80)
#pl.show()

#print "Eigenvalues", evals
#print "Eigenvectors", evecs
#print "first Eigenvector..."
#print evecs[0]

fig=pl.figure()

#pl.title("DoS by TightBinding")
#pl.subplot(311) #3 subplots stacked on top of one another

#Plot Eigenvalues with filled-in Eigenvectors^2 / electron probabilities
#pl.subplot(311)
#for j in [0,n/2]: #range(0,5): #Number of eigenvalues plotted (electron wavefns)
#    psi=evecs[:,j]*evecs[:,j]
#    pl.fill_between(range(n),0,psi, facecolor=colours[j%8])
#    pl.plot(range(n),evecs[:,j],color=colours[j%8])
#pl.ylabel("Occupation")
#pl.ylim((3.8,5.0))
#pl.yticks(fontsize=9)
#pl.xticks(visible=False)

#Plot cumulative eigenvectors / probability density
#pl.subplot(312)
#for j in [0,1,2,n/4, n/2]: #range(0,5): #Number of eigenvalues plotted (electron wavefns)
#    psi=evecs[:,j]*evecs[:,j]    # expectation value
#    pl.fill_between(range(n),0,sorted(psi,reverse=True), facecolor=colours[j%8]) #can't see this anymore on large plots...
#    psi=sorted(psi,reverse=True) # expectation value, ranked in order (largest first)

#    psi_sum=[0.0]
#    for i in range(len(psi)): # should be a nicer way to do this with functional programming!
#        psi_sum.append(psi_sum[-1]+psi[i])
    
    #    pl.plot(psi_sum, color=colours[j%8])
#    pl.plot(y=0.95) # 2 sigma confidence interval? # TODO: why doesn't this work?
#pl.ylabel("Cumulative Density")
#pl.ylim((3.8,5.0))
#pl.yticks(fontsize=9)
#pl.xticks(visible=False)


#Plot DoS
pl.hist(evals,100,histtype='stepfilled',color='b')
#pl.ylim(0,80)
pl.xlabel("Energy (eV)")
pl.ylabel("DOS")
pl.yticks(fontsize=9)
pl.xlim(-3.85,-3.45)

pl.tight_layout(pad=0.3)

pl.show() #Displays plots!

print "Lowest Eigenvalue:\n", evals[0]

fp=open('eigenvector_balls_pymol.py','w')
fp.write("from pymol.cgo import *    # get constants \nfrom pymol import cmd \n")

psi=evecs[:,0]*evecs[:,0] # scale everything relative to max density on first eigenvector

for ei,colour in zip( [0,5,10,50] , [(0,0,1),(0,1,0),(1,1,0),(1,0,0)]):
    psi=evecs[:,ei]*evecs[:,ei]
    maxpsi=max(psi)
    
    fp.write("obj = [\n")

    for i in reversed(np.argsort(psi)): #magic list of sorted array indices
        #       print locations[i]
        #       print psi[i]
        weight=float(psi[i])/maxpsi #on interval [0,1]
        fp.write("ALPHA, %f,\n" %(weight))
        weight=1
        fp.write("COLOR, %f, %f, %f,\n" %(colour[0]*weight , colour[1]*weight, colour[2]*weight))
        fp.write("SPHERE, %f, %f, %f, 5.0,\n" %(filein_pos[i][1],filein_pos[i][2],filein_pos[i][3]))
    
    fp.write(" END \n]\ncmd.load_cgo(obj,'EV_%d') \n" %(ei))



#print r

#print "Effective size of polaron is from", np.min(r)
#print "to", np.max(r)

#print "Saving figures...(one moment please)"
now=datetime.datetime.now().strftime("%Y-%m-%d-%Hh%Mm") #String of standardised year-leading time
pl.annotate("%s"%now,xy=(0.75,0.02),xycoords='figure fraction') #Date Stamp in corner

fig1.savefig("%s-polaron.pdf"%now)
fig.savefig("%s-DBTW.pdf"%now) #Save figures as both PDF and easy viewing PNG (perfect for talks)
#fig.savefig("%s-DBTW.png"%now)
#fig.savefig("%s-LongSnakeMoan.ps"%now)    # TODO: check with latest python scripts to see best way to export these for future inclusion in Latex etc.

