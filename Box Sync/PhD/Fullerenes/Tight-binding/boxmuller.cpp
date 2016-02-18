// cg-pcbm-j.c
//  Generate Coarse grain Js - calculate isotropic Js from Coarse Grain PCBM or
//  C60 molecular dynamics simulations
// File originally called 'gen_edges.c'; Jarvist Moore Frost 2013

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

// for size of data structure. Nb: check your ulimits for heap size if
// segfault!
#define MAXPOINTS 100000


double uniformRandom()
{
    return ( (double)(rand())+1.)/((double)(RAND_MAX)+1.);
}


double gauss(double dx)
{
    double noise;
    double u1=uniformRandom();
    double u2=uniformRandom();
    
    noise=dx*sqrt(-2.*log(u1))*cos((u2)*2*M_PI);
    
    return noise;
}


int main()
{
    int i;
    double rands;
    
    for (i=0; i<100; i++){
        rands = gauss(0.1);
        cout<<rands<<endl;
    }
    
    return 0;

    
}
