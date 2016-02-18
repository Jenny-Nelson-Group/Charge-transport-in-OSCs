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

double uniformRandom()
{
    return ( (double)(rand())+1.)/((double)(RAND_MAX)+1.);
}


double gauss(double x, double dx)
{
    double noise;
    double u1=uniformRandom();
    double u2=uniformRandom();
    
    noise=sqrt(dx*-2.*log(u1))*cos(u2)*2*M_PI;
    
    x=x+noise;
    
    return x;
}

int main()
{
    double coords[MAXPOINTS][3];
    double d,mind,r[3],s;
    int n,i,j;
    
    double lambda=0.6;     // DRAGONS!
    double prefactor=10.0; //check this!
    
    n=4;
    mind=30.0;
    
    if (n>MAXPOINTS)
    fprintf(stderr,"ERROR: Too many points! Will probably now crash...");
    
    for (i=0;i<n;i++) // read tuple coordinates in
    scanf("%lf %lf %lf", & coords[i][0], & coords[i][1], & coords[i][2]);
    
//int k;
//double box_r = 1.0/box;
//   double box2 = 0.5*box;
//   double box2_r = 1.0/box2;
    
    for (i=0;i<n;i++) // N^2 search over particle pairs (!!!) - this is why we are using C
    for (j=0;j<i;j++) // Nb: limit to i to just do lower diagonal of n*n
    if (i!=j) // no infinities if self
    {
        //r[0]=coords[i][0]-coords[j][0];
        //r[1]=coords[i][1]-coords[j][1];
        //r[2]=coords[i][2]-coords[j][2];
        
        r[0]=coords[i][0]-coords[j][0];           // Second line is PBCs
//    r[0]-=round(r[0]/n)*n;
        
        r[1]=coords[i][1]-coords[j][1];
     //   r[1]-=round(r[1]/n)*n;
        
        r[2]=coords[i][2]-coords[j][2];
       // r[2]-=round(r[2]/n)*n;

        
        d=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
        
        
        s=-6.0;
        s=gauss(s,1);
        
        //  printf("%d %d \t%f %f %f\t%f\n",i,j,coords[i][0],coords[i][1],coords[i][2],d);
        if (d<=mind)
        printf("%d\t%d\t%lf\t%lf\n",i,j,prefactor*exp(-d*lambda),s);
        //increasing factor makes J tail off 'sharper'
    }
    
}
