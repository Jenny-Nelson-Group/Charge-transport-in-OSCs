// cg-pcbm-j.c
//  Generate Coarse grain Js - calculate isotropic Js from Coarse Grain PCBM or
//  C60 molecular dynamics simulations
// File originally called 'gen_edges.c'; Jarvist Moore Frost 2013

#include <stdio.h>
#include <math.h>

// for size of data structure. Nb: check your ulimits for heap size if
// segfault!
#define MAXPOINTS 100000

int main()
{
    double coords[MAXPOINTS][3];
    double d,mind,r[3];
    int n,i,j;
    
    double lambda=6;     // DRAGONS!
    double prefactor=10.0; //check this!
    
    n=1000;
    mind=1.0;
    
    
    if (n>MAXPOINTS)
    fprintf(stderr,"ERROR: Too many points! Will probably now crash...");
    
    for (i=0;i<n;i++) // read tuple coordinates in
    scanf("%lf %lf %lf",& coords[i][0], & coords[i][1], & coords[i][2]);
    
    for(i=0;i<n;i++){
        coords[i][0]-=floor(coords[i][0]/n)*n;    //PBCs
        coords[i][1]-=floor(coords[i][1]/n)*n;
        coords[i][2]-=floor(coords[i][2]/n)*n;
    }
    
    for (i=0;i<n;i++) // N^2 search over particle pairs (!!!) - this is why we are using C
    for (j=0;j<i;j++) // Nb: limit to i to just do lower diagonal of n*n
    if (j!=i) // no infinities if self
    {
        r[0]=coords[i][0]-coords[j][0];           // Second line is PBCs
        r[0]-=rint(r[0]/n)*n;
        
        r[1]=coords[i][1]-coords[j][1];
        r[1]-=rint(r[1]/n)*n;
        
        r[2]=coords[i][2]-coords[j][2];
        r[2]-=rint(r[2]/n)*n;
        
        d= sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        
        //  printf("%d %d \t%f %f %f\t%f\n",i,j,coords[i][0],coords[i][1],coords[i][2],d);
        if (d<=mind)
        printf("%d\t%d\t%lf\n",i,j,prefactor*exp(-d*lambda));
        //increasing factor makes J tail off 'sharper'
    }
    
}
