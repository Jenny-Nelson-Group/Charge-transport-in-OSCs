// cg-pcbm-j.c
//  Generate Coarse grain Js - calculate isotropic Js from Coarse Grain PCBM or
//  C60 molecular dynamics simulations
// File originally called 'gen_edges.c'; Jarvist Moore Frost 2013

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// for size of data structure. Nb: check your ulimits for heap size if
// segfault!
#define MAXPOINTS 10000
#define TWO_PI 6.2831853071795864769252866
#define n 5


double uniformRandom()
{
    return ( (double)(rand())+1.)/((double)(RAND_MAX)+1.);
}


double gauss(double x, double dx)            // Function to add gaussian disorder
{
    double noise;
    double u1=uniformRandom();
    double u2=uniformRandom();
    
    noise=sqrt(dx*-2.*log(u1))*cos(u2)*TWO_PI;
    
    x+=noise;
    
    return x;
}


int main()
{
    double coords[MAXPOINTS][3];
    double d[n],d_dis[n],mind,r[3];
    int i,j;
    
    double lambda=6.;     // DRAGONS!
    double prefactor=10.0; //check this!
    
    mind=1.;
    
    if (n>MAXPOINTS)
    fprintf(stderr,"ERROR: Too many points! Will probably now crash...");
    
    for (i=0;i<n;i++) // read tuple coordinates in
    scanf("%lf %lf %lf", & coords[i][0], & coords[i][1], & coords[i][2]);
    
    
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
        r[0]-=round(r[0]/n)*n;
        
        r[1]=coords[i][1]-coords[j][1];
        r[1]-=round(r[1]/n)*n;
        
        r[2]=coords[i][2]-coords[j][2];
        r[2]-=round(r[2]/n)*n;
        
        d[i]=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        
        d_dis[i]=gauss(d[i],0);     // Adding disorder with variance as argument
        
        
        //  printf("%d %d \t%f %f %f\t%f\n",i,j,coords[i][0],coords[i][1],coords[i][2],d);
        if (d[i]<=mind)
        //printf("%d\t%d\t%lf\t%lf\n",i,j,prefactor*exp(-d[i]*lambda),prefactor*exp(-d_dis[i]*lambda));
        printf("%d\t%d\t%lf\n",i,j,prefactor*exp(-d_dis[i]*lambda));
        //increasing factor makes J tail off 'sharper'
    }
    
}
