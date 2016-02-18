// cg-pcbm-j.c
//  Generate Coarse grain Js - calculate isotropic Js from Coarse Grain PCBM or
//  C60 molecular dynamics simulations
// File originally called 'gen_edges.c'; Jarvist Moore Frost 2013

#include <stdio.h>
#include <math.h>

// for size of data structure. Nb: check your ulimits for heap size if
// segfault!
#define MAXPOINTS 10000

int main()
{
    double coords[MAXPOINTS];
    double d,mind,r[1];
    int n,i,j;
    
    double lambda=0.6;     // DRAGONS!
    double prefactor=10.0; //check this!
    
    n=1000;
    mind=30;
    
    if (n>MAXPOINTS)
    fprintf(stderr,"ERROR: Too many points! Will probably now crash...");
    
    for (i=0;i<n;i++) // read tuple coordinates in
    scanf("%lf", & coords[i]);

    
    for (i=0;i<n;i++) // N^2 search over particle pairs (!!!) - this is why we are using C
    for (j=0;j<i;j++) // Nb: limit to i to just do lower diagonal of n*n
    if (j!=i) // no infinities if self
    {
        d=coords[i]-coords[j];
        
        //  printf("%d %d \t%f %f %f\t%f\n",i,j,coords[i][0],coords[i][1],coords[i][2],d);
        if (d<=mind)
        printf("%d\t%d\t%lf\n",i,j,prefactor*exp(-d*lambda));
        //increasing factor makes J tail off 'sharper'
    }
    
}
