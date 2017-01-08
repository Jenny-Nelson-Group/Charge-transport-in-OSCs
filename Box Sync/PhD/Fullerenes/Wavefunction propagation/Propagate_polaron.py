#!/usr/bin/env python

# Read in polaron wavefunction and propagate with TD Schroedinger equation

import numpy as np
import matplotlib.pyplot as pl
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import sys
import matplotlib.animation as animation
from sys import exit
from IPython import embed
import pickle
import warnings
warnings.simplefilter("ignore", np.ComplexWarning)


#---------------------- Function to calculate occupation at a given time. Reads in time, evals and evecs of unperturbed (TB) Hamiltonian, and evecs of the polaron wavefunction. Reads out bins of increasing radius along with the total occupation within that radius-----------------------------

def calcocc(T,Evals,Evecs_init,Pvecs_polaron):
    
    psi1=Pvecs_polaron      # Polaron wavefunction to propagate
    centre=np.zeros(3)
    psi1_occ=np.zeros(n)
    r=np.zeros(n)
    dist_from_centre=np.zeros((n,3))

    Evecs=np.array(Evecs_init,dtype=complex)
   
 
    for i in range(0,n):
        psi1_occ[i]=psi1[i]*np.conj(psi1[i]) # Calculate occupations of each site with polaron
    
    max_index=np.argmax(psi1_occ)     # Find which site is maximally occupied to use as centre
    centre[0]=locations[max_index,0]
    centre[1]=locations[max_index,1]
    centre[2]=locations[max_index,2]

#print "centre: ", centre

#Positions of all molecules from max occupied one with PBCs

    for j in range(0,n):
        dist_from_centre[j,:]=locations[j,:]-centre

    print dist_from_centre
	
    dist_from_centre=dist_from_centre/cell
    dist_from_centre[dist_from_centre>0.5]-=1.0
    dist_from_centre=dist_from_centre*cell
	
    for j in range(0,n):
        r[j]=np.linalg.norm(dist_from_centre[j,:])    

    idx=np.argsort(r)      #Sort distances
    sort_r=r[idx]
	
    print r

    bins=np.linspace(np.amin(sort_r),np.amax(sort_r),nbins)     # Bins for different radii 

    inds=np.digitize(sort_r,bins)                  # Indices for each r depending on bin      
    
    hbar=6.582*10**-16 # eV.s 
    z=1j               # imaginary unit         
    
    psi_t=np.zeros(n,dtype=complex)
    psi_t_occ=np.zeros(n,dtype=complex)
    psi_bin=np.zeros(nbins,dtype=complex)
    
    for i in range(0,n):                  # Find coefficient a and wavefunction after time T 
        a=np.inner(Evecs[:,i],psi1)		  # with TD Schrodinger equation
        psi_t+=a*np.exp(-z*Evals[i]*dt*T/hbar)*Evecs[:,i]
    
    psi_t_norm=np.linalg.norm(psi_t)      # Normalise
    psi_t=psi_t/psi_t_norm

    for j in range(0,n):
        psi_t_occ[j]=psi_t[j]*np.conj(psi_t[j])   # Occupations after time T

    psi_t_occ=psi_t_occ.astype(np.float_)

    psi_occ=psi_t_occ[idx]                        # Sort occupations by radius from centre

    for i in range(0,n):                          # Save occupations in bins depending            
        psi_bin[inds[i]-1]+=psi_occ[i]            # on radius from centre  
   

    psi_bin=psi_bin.astype(np.float_)

   
    return psi_t, bins, psi_bin


#------ Plot snapshots of occupation vs distance from centre at various times on a 3D plot -------

def plotpropagation():
    
    tsteps=10              #No. of plots
    
    Z=np.zeros((nbins,tsteps))
    X_all=np.zeros((nbins,n))
    X=np.zeros(nbins)
    T=np.linspace(0,tsteps,tsteps)
    
    for t in range(0,tsteps):
        for i in range(0,n):
            psi_t,X_all[:,i],Z[:,t]=calcocc(t,evals,evecs,psis[:,i])   # Define the variables (X is radius bins
    													# Z[:,t] is occupation at time t)
	for i in range(0,nbins):
		X[i]=np.mean(X_all[i,:])

    print X
    
    fig=pl.figure()           #Plot time evolved wavefunction for different dts
    
    ax = fig.add_subplot(111, projection='3d')
    
    verts = []
    for i in xrange(T.shape[0]):                        # Define lines to plot on 3D plot  
        verts.append(zip(X[0:100],Z[:,i]))
    
    poly = PolyCollection(verts, facecolors=(1,1,1,1), edgecolors=(0,0,1,1))
    ax.add_collection3d(poly, zs=T, zdir='y')
    ax.set_xlim3d(np.min(X),100)
    ax.set_ylim3d(np.min(T), np.max(T))
    ax.set_zlim3d(np.min(Z), np.max(Z))
    ax.set_xlabel('Radius from centre (A)')
    ax.set_ylabel('Time (fs)')
    ax.set_zlabel('Occupation')
    ax.view_init(20,120)
    
    fig.savefig("WF_propagation_%s_allpolarons.pdf"%(coordfile))


#----------------- Animate occupation versus distance from centre over time ----------------------

def Animate():                              
    
	sort_r_all=np.zeros((nbins,n))	
	sort_r=np.zeros(nbins)
	fig=pl.figure()
    
	xlim=(0,60)
	ylim=(0,1)
    
	ax = pl.axes(xlim=xlim,ylim=ylim)
	psi_x_line, = ax.plot([],[],c='r')
    
	def init():
		psi_x_line.set_data([], [])       # Initial plot 
		return psi_x_line,
	print "Calculating parameters..."
    
	def animate(i):
		for j in range(0,n):
			psi_t,sort_r_all[:,j],psi_occ=calcocc(i,evals,evecs,psis[:,j])  # Calculate occupation vs radius for
		for j in range(0,n):   
			sort_r[i]=np.mean(sort_r_all[i,:])	    	 	
		
		x=sort_r											 # different times i 
 		y=psi_occ
		psi_x_line.set_data(x,y)   
 		return psi_x_line,
	print "Animating..."
    
	anim=animation.FuncAnimation(fig,animate,init_func=init,frames=1000,interval=10,blit=False)
    
	pl.show()
	print "Saving..."
	anim.save('wavefunction_propagation_allpolarons_%s.mp4'%(coordfile),fps=30,writer='ffmpeg',extra_args=['-vcodec', 'libx264'])

# Output r and root_t for finding diffusion constant

def r_vs_root_t(Evecs_polaron):
    psi_t,BINS,PSI_BIN=calcocc(Time,evals,evecs,Evecs_polaron)

    cum_occ=np.cumsum(PSI_BIN)     # Cumulative sum for occupation to find radius containing 99% of charge

    for i in range(0,n):
        if cum_occ[i]>0.99:break
    r_c=BINS[i]

    root_t=np.sqrt(Time*dt)

    return root_t, r_c



#---------------------- Save stills of polaron to visualise in pymol ------------------------------

def polaronvisualise(Psi_t):
    fp=open('Pymol_%s_%s.py'%(coordfile,Time),'w')
    fp.write("from pymol.cgo import *    # get constants \nfrom pymol import cmd \n")
    
    psi=Psi_t*Psi_t # scale everything relative to max density on first eigenvector
   
    maxpsi=max(psi)
    colour= (0,0,1)   
    fp.write("obj = [\n")
        
    psisum=0.0
    for i in reversed(np.argsort(psi)): #magic list of sorted array indices
        weight=float(psi[i])/maxpsi #on interval [0,1]
        fp.write("ALPHA, %f,\n" %(weight))
        weight=1
        fp.write("COLOR, %f, %f, %f,\n" %(colour[0],colour[1],colour[2]))
        fp.write("SPHERE, %f, %f, %f, 5.0,\n" %(locations[i][0],locations[i][1],locations[i][2]))
        
    fp.write(" END \n]\ncmd.load_cgo(obj,'%s_%s') \ncmd.bg_color(color='white') \ncmd.zoom(complete=1) \n"%(coordfile,Time))



# ----------------------------------- Read in files -----------------------------------------

evals=np.loadtxt(sys.argv[1])
evecs=np.loadtxt(sys.argv[2])
coordfile=sys.argv[3]
n=int(sys.argv[4])
Time=float(sys.argv[5])
adduct=sys.argv[6]
cell=[float(sys.argv[7]),float(sys.argv[8]),float(sys.argv[9])]

dt=1*10**-15           # Timestep for plots
nbins=100			   # Number of radial bins

locations=np.loadtxt(coordfile)

psis=np.zeros((n,n))

for i in np.arange(0,n,50):
    num_string="00%i"%(i)
    psi=np.loadtxt("Evecs_%s_%s_polaron.dat"%(adduct,num_string[-3:]))
    psi=psi.reshape((n,50))
    psis[:,i:i+50]=psi


#----------------------------------- Call functions -------------------------------------------

#Animate()
#plotpropagation()

calcocc(Time,evals,evecs,psis[:,0])

rs=np.zeros(n)


# Including all polarons
#for i in range(0,n):
#   sqrt_t,rs[i]=r_vs_root_t(psis[:,i])

#R_C=np.mean(rs)

# Including only polarons not at the edges

bulk_indices=[]

#for i in range(0,n):
#	if 30<locations[i,0]<(cell[0]-30) and 30<locations[i,1]<(cell[1]-30) and 30<locations[i,2]<(cell[2]-30):
#		bulk_indices=np.append(bulk_indices,i)

#bulk_indices = [int(x) for x in bulk_indices]           # Set elements as integers

#m=len(bulk_indices)
#rs_noedge=np.zeros(m)

#psis_noedge=psis[:,bulk_indices]

#for i in range(0,m):
#   sqrt_t,rs_noedge[i]=r_vs_root_t(psis_noedge[:,i])

#R_C=np.mean(rs_noedge)

#print sqrt_t,R_C



#polaronvisualise(psi_t)











