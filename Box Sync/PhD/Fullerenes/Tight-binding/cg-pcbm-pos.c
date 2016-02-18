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
    int n,i;
    
    n=1000;
    
    if (n>MAXPOINTS)
        fprintf(stderr,"ERROR: Too many points! Will probably now crash...");
    
    for (i=0;i<n;i++) {// read tuple coordinates in
    
    
        scanf("%lf %lf %lf",& coords[i][0], & coords[i][1], & coords[i][2]);
    
        printf("%d\t%f\t%f\t%f\n",i,coords[i][0],coords[i][1],coords[i][2]);
    }

    
}
