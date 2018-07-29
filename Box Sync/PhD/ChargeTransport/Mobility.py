#!/usr/bin/env python

""" Calculate the mobility given a crystal structure and the electronic couplings (Js) from the central unit cell to all others in a 3x3x3 supercell. The J file should be in the format $i $j J, where $i and $j are the indices for each molecule, and this code is set up the take advantage of the symmetry in a crystal and duplicate the Js between additional unit cells. The code Transform_molecules.py and can be used to make up the unit cells, given symmetry operations for the crystal structure.
    The code works by using the populated J matrix to find the corresponding rates in the prescence of an electric field (although field strength can be set to zero and the mobility calculated via the Einstein relation).
    A master equation for rates is then set up and solved: AP=0, where A is the rate matrix containing the sum of all rates away from a molecule on the diagonal and the rate to each other molecule on the off diagonals. This is used to find the velocity, which can then be used to find the mobility via the field.
    The diffusion coefficient is calculated with a model based on a random walk, which is converted to a mobility with the Einstein relation. """


import numpy as np
import matplotlib.pyplot as pl
import scipy
from scipy import linalg, matrix
from scipy.sparse.linalg import svds
import sys
from sympy import Matrix


LAMBDA=0.3                     # eV
HBAR=6.582*10**-16             # eV.s
e=1                            # Charge on electron in eV/V
KB=8.617*10**-5                # in eV/K
T=300                          # in K


def null(X,eps=1e-8):
    """ Define null space function for solving AP=0 """
    Solution=0
    u, s, vh = scipy.linalg.svd(X)     # Single value decomposition
    null_mask = (np.absolute(s) <= eps)
    null_space = scipy.compress(null_mask, vh, axis=0)  # Find corresponding vectors
    n,m=np.shape(null_space)
    for i in range(0,n):           # Choose only solutions that can be probabilities
        positive=np.greater_equal(null_space[i,:],np.zeros(m))
        negative=np.less_equal(null_space[i,:],np.zeros(m))
        if (np.all(positive)==True or np.all(negative)==True):        # Ensure sum(Pi)=1
            Solution=np.absolute(null_space[i,:])
            Solution=Solution/np.sum(Solution)
    return scipy.transpose(Solution)


def find_rows(a,b):
    """ Find matching rows for populating full J matrix using crystal symmetry """
    
    row_match = np.array(np.all((np.isclose(a[:,None,:],b[None,:,:],rtol=1e-3,atol=1e-8)),
                                axis=-1).nonzero()).T.tolist()
                                
    return row_match


def fill_Js(J,distancematrix,N,M):
    """ Find sites in 1st unit cell (by finding rows of matrix containing more than a critical value of non-zeros) """
    
    unit_cell_sites=[]
    
    for i in range(0,N):
        if len(J[np.nonzero(J[i,:])])>M+1:
            unit_cell_sites=np.append(unit_cell_sites,[i])

    unit_cell_sites = [int(x) for x in unit_cell_sites]      # Set elements as integers

    unit_cell_vecs=np.zeros((N,N,3))
    rest_distance_matrix=np.zeros((N,N,3))

    unit_cell_vecs[0:M,:,:]=distancematrix[0:M,:,:]
    rest_distance_matrix[M:N,:,:]=distancematrix[M:N,:,:]

    for i in range(0,M):
        for j in range(M,N):
            J_match = find_rows(unit_cell_vecs[i,:,:],rest_distance_matrix[j,:,:])
            if unit_cell_vecs[i,J_match[0][0],:].all()!=0:
                for match in range(0,len(J_match)):
                    if i!=match:
                        J[j,J_match[match][1]]=J[i,J_match[match][0]]
                        J[J_match[match][1],j]=J[J_match[match][0],i]

    print len(J[np.nonzero(J)])

    return J


def Marcus_rates(J,field,distancematrix,N):
    """ Calulate matrix of rates with Marcus theory """
    
    A=np.zeros((N,N))
    
    for i in range(0,N):
        for j in range(0,N):
            if (i!=j):
                deltaE=0
                d=distancematrix[i,j,:]
                Field_d=np.dot(field,d)
                A[i,j]=((J[i,j]*J[i,j]))*((np.pi/(LAMBDA*KB*T))**0.5)*np.exp(-(((deltaE-Field_d)+LAMBDA)**2)/(4*LAMBDA*KB*T))
    return A


def Master_eq(A,N):
    """ Solve Master equation to find steady state probablities """
    
    for i in range(0,N):
        A[i,i]=-np.sum(A[:,i])
    
    P_all=null(A)

    return P_all


def mob(A,P_all,field_unit,F_mag,distancematrix,N):
    """ Find mobility by calculating velocity. Can only be used for field strength != 0 """
    
    V=np.zeros(3)       # Find velocity vector
    
    for i in range(0,N):
        for j in range(0,N):
            if i!=j:
                V+=distancematrix[i,j,:]*A[i,j]*P_all[i]


    V_F=np.dot(V,field_unit)       # V in direction of field

    mu = (V_F/F_mag)/HBAR                # Mobility in cm^2/Vs

    print "Mobility= ", mu, " cm^2/Vs"
    
    return mu


def mob_einstein(A,P_all,field,distancematrix,N):
    """ Find mobility from Einstein relation. Can be used even for zero field. """
    
    # Calculate D. Initialise.
    D = 0
    
    for i in range(0,N):
        for j in range(0,N):
            D += 0.5*P_all[i]*A[j,i]*(np.dot(distancematrix[j,i,:],field))**2

    mu_ein = ((e*D)/(KB*T))/HBAR

    print "Mobility from Einstein relation= ", mu_ein, " cm^2/Vs"
    
    return mu_ein


def change_field(Js_all,F_mag,distancematrix,N,mob_type):
    """ Apply field at 10 degree angles in each of the 3 orthogonal planes """
    
    mobilities_ab = []
    mobilities_bc = []
    mobilities_ac = []
    
    for angle in np.linspace(0,2*np.pi,36):
        Field_ab=[F_mag*np.cos(angle),F_mag*np.sin(angle),0]
        Field_bc=[0,F_mag*np.cos(angle),F_mag*np.sin(angle)]
        Field_ac=[F_mag*np.cos(angle),0,F_mag*np.sin(angle)]
        
        A_ab = Marcus_rates(Js_all,Field_ab,distancematrix,N)
        A_bc = Marcus_rates(Js_all,Field_bc,distancematrix,N)
        A_ac = Marcus_rates(Js_all,Field_ac,distancematrix,N)
        
        P_ab = Master_eq(A_ab,N)
        P_bc = Master_eq(A_bc,N)
        P_ac = Master_eq(A_ac,N)
        
        if mob_type=='mob':
            mob_ab = mob(A_ab,P_ab,Field_ab/np.linalg.norm(Field_ab),F_mag,distancematrix,N)
            mob_bc = mob(A_bc,P_bc,Field_bc/np.linalg.norm(Field_bc),F_mag,distancematrix,N)
            mob_ac = mob(A_ac,P_ac,Field_ac/np.linalg.norm(Field_ac),F_mag,distancematrix,N)
        
        if mob_type=='mob_einstein':
            mob_ab = mob_einstein(A_ab,P_ab,Field_ab/np.linalg.norm(Field_ab),distancematrix,N)
            mob_bc = mob_einstein(A_bc,P_bc,Field_bc/np.linalg.norm(Field_bc),distancematrix,N)
            mob_ac = mob_einstein(A_ac,P_ac,Field_ac/np.linalg.norm(Field_ac),distancematrix,N)
        
        mobilities_ab = np.append(mobilities_ab,mob_ab)
        mobilities_bc = np.append(mobilities_bc,mob_bc)
        mobilities_ac = np.append(mobilities_ac,mob_ac)

    return np.array(mobilities_ab),np.array(mobilities_bc),np.array(mobilities_ac)


def plot_mobility(mobilities_ab,mobilities_bc,mobilities_ac,J_cutoff,filename):
    """ Plot the mobility on a polar plot """
    
    np.savetxt("Mob_%s_%s_%s.txt"%('ab',filename,J_cutoff),mobilities_ab)
    np.savetxt("Mob_%s_%s_%s.txt"%('bc',filename,J_cutoff),mobilities_bc)
    np.savetxt("Mob_%s_%s_%s.txt"%('ac',filename,J_cutoff),mobilities_ac)
    
    max_mob = np.max([mobilities_ab,mobilities_bc,mobilities_ac])
    min_mob = np.min([mobilities_ab,mobilities_bc,mobilities_ac])
    av_mob = np.mean([mobilities_ab,mobilities_bc,mobilities_ac])
    
    print "Max mobility: ", max_mob
    print "Min mobility: ", min_mob
    print "Average mobility: ", av_mob
    
    fig1=pl.figure()
    
    ax=pl.subplot(111, polar=True)
    ax.plot(np.linspace(0,2*np.pi,36),mobilities_ab,'o',label='Mobility (cm$^2$V$^{-1}$s$^{-1}$)')
    ax.grid(True)
    ax.set_rmax(max_mob)
    ax.legend(loc='upper left', bbox_to_anchor=(-0.2,1.1))
    
    #pl.show()
    #fig1.savefig("%s_mobility_ab_plane.png"%(filename))
    
    fig2=pl.figure()
    
    ax=pl.subplot(111, polar=True)
    ax.plot(np.linspace(0,2*np.pi,36),mobilities_bc,'o',label='Mobility (cm$^2$V$^{-1}$s$^{-1}$)')
    ax.grid(True)
    ax.set_rmax(max_mob)
    ax.legend(loc='upper left', bbox_to_anchor=(-0.2,1.1))
    
    #pl.show()
    #fig2.savefig("%s_mobility_bc_plane.png"%(filename))
    
    fig3=pl.figure()
    
    ax=pl.subplot(111, polar=True)
    ax.plot(np.linspace(0,2*np.pi,36),mobilities_ac,'o',label='Mobility (cm$^2$V$^{-1}$s$^{-1}$)')
    ax.grid(True)
    ax.set_rmax(max_mob)
    ax.legend(loc='upper left', bbox_to_anchor=(-0.2,1.1))
    
    #pl.show()
    #fig3.savefig("%s_mobility_ac_plane.png"%(filename))
    
    return 0


def main():
    """ Load data and run functions """
    
    coordfile = np.loadtxt(sys.argv[1]) # read in coordinate file (in Angstroms)
    Jif = np.loadtxt(sys.argv[2])          # Read in J file (with columns $i $j J)
    F_mag = float(sys.argv[3])          # Magnitude of field vector (in V/cm)
    J_cutoff=float(sys.argv[4])       # Define a cutoff for Js to ignore
    filename = sys.argv[5]
    
    N = len(coordfile)                  # Number of molecules
    
    M = N/27                         # M = number of molecules in unit cell (for 3x3x3 supercell)
    
    t = 2

    Jif=Jif[Jif[:,2]>J_cutoff]

    Size = len(Jif[:,2])

    orderedJs = np.argsort(Jif[:,t])

    print "Top Js: ", Jif[orderedJs[Size-10:Size],2], "at", Jif[orderedJs[Size-10:Size],0], Jif[orderedJs[Size-10:Size],1]

    print "Max J: ", np.max(Jif[:,t])

    Js = np.zeros((N,N))              # Set up Js in matrix

    for i in range(0,Size):
        Js[int(Jif[i,0]),int(Jif[i,1])]=Jif[i,t]
        Js[int(Jif[i,1]),int(Jif[i,0])]=Jif[i,t]
    
    coordfile = coordfile*10**-8     # Convert from Angstroms to cm
    
    distancematrix = coordfile[:,None,...]-coordfile[None,...]
    
    J_all = fill_Js(Js,distancematrix,N,M)
    
    mobilities_ab,mobilities_bc,mobilities_ac = change_field(J_all,F_mag,distancematrix,N,'mob')
    
    plot_mobility(mobilities_ab,mobilities_bc,mobilities_ac,J_cutoff,filename)


main()
