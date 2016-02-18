// cg-pcbm-j.c
//  Generate Coarse grain Js - calculate isotropic Js from Coarse Grain PCBM or
//  C60 molecular dynamics simulations
// File originally called 'gen_edges.c'; Jarvist Moore Frost 2013

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// for size of data structure. Nb: check your ulimits for heap size if
// segfault!
#define MAXPOINTS 100000



double uniformRandom()
{
    return ( (double)(rand())+1.)/((double)(RAND_MAX)+1.);
}


double gauss(double x, double dx)
{
    double noise;
    double u1=uniformRandom();
    double u2=uniformRandom();
    
    noise=dx*sqrt(-2.*log(u1))*cos((u2)*2*M_PI);
    
    x+=noise;
    
    return x;
}


int main()
{
    double coords[MAXPOINTS][3];
    double d,s,s_dis,mind,r[3];
    double box, box_r, box2, box2_r;
    int n,i,j,k,p;
    double x,y,z;
    
    double lambda=0.6;     // DRAGONS!
    double prefactor=10; //check this!
    
    scanf("%lf %d %lf", &box,&n,&mind); // First line of STDIN is 'box side length' 'number of points' 'min distance to consider'
    
    if (n>MAXPOINTS)
        fprintf(stderr,"ERROR: Too many points! Will probably now crash...");
    
    //PBCs
    
    box_r=1.0/box;
    box2=0.5*box;
    box2_r=1.0/box2;
    
    for (i=0;i<n;i++) // read tuple coordinates in
        scanf("%lf %lf %lf",& coords[i][0], & coords[i][1], & coords[i][2]);
    
    
    for (i=0;i<n;i++) // N^2 search over particle pairs (!!!) - this is why we are using C
        for (j=0;j<i;j++) // Nb: limit to i to just do lower diagonal of n*n
            if (j!=i) // no infinities if self
            {
                for (k=0;k<3;k++)
                r[k]=coords[i][k]-coords[j][k];
                p=(int)r[k]*box2_r;
                r[k]-=p*box;
                
                x=r[0];
                y=r[1];
                z=r[2];
                
                d=sqrt(x*x+y*y+z*z);
                
                s=-3.7;
                s_dis=gauss(s,0);
                
                if (d<=mind)
                    printf("%d\t%d\t%lf\t%lf\n",i,j,prefactor*exp(-d*lambda),s_dis); //increasing factor makes J tail off 'sharper'
            }
    
}
