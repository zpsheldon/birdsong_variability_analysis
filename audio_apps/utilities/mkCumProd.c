/*=================================================================
 *
 * mkCumProd.c	
 *	        cumulative outerproduct for dynamic time warping
 *          paths.
 *
 * The calling syntax is:
 *
 *		D = mkCumProd(d,pwt,mdim,ndim,mlim)
 *      d: raw outerproduct
 *      pwt: path weightings
 *      mdim/ndim: d rows/cols 
 *      mlim: limit defining Sakoe-Chiba band (Rabiner & Juang 1993)  
 *
 *      values outside the band are 0 by default    
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2000 The MathWorks, Inc.
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"

void mkCumProd(
		   double	*d,
           double   *dCum,
           double  *pwt,
           double mdim,
           double ndim,
           double mlim
		   )
{

    double	val1,val2,val3;
    int n,m,lm1,lm2,p11,p12,p23,p21,p22,p32;

    for (n = 3; n <= ndim; n++) {
        lm1 = n-mlim;
        lm2 = mlim+n-1;
        lm1 = lm1>3? lm1:3;
        lm2 = lm2<mdim? lm2:mdim;

        for (m = lm1; m <= lm2; m++) {
            
            p11 = n+(m-1)*ndim-1;
            
            p21 = n-1+(m-1)*ndim-1;
            p32 = n-2+(m-2)*ndim-1;
            
            p12 = n+(m-2)*ndim-1;
            p23 = n-1+(m-3)*ndim-1;
            
            p22 = n-1+(m-2)*ndim-1;
            
            val1 = *(pwt) * (.5 * *(d+p11) + .25 * *(d+p21) + .25 * *(d+p22)) + *(dCum+p32);      
            val2 = *(pwt+1) * (.5 * *(d+p11) + .25 * *(d+p12) + .25 * *(d+p22)) + *(dCum+p23);
            val3 = *(pwt+2) * *(d+p11) + *(dCum+p22);
            
            *(dCum+p11) = val1>val2? (val1>val3? val1:val3) : (val2>val3? val2:val3);
            
        }
        
    }
    
}

void mexFunction( int nlhs, mxArray*plhs[], 
		  int nrhs, const mxArray*prhs[])
     
{ 
    double *dCum,*d;
    double *pwt;
    int mdim,ndim,mlim;
    
    d = mxGetPr(prhs[0]);
    pwt = mxGetPr(prhs[1]);
    
    mdim = mxGetScalar(prhs[2]);
    ndim = mxGetScalar(prhs[3]);
    mlim = mxGetScalar(prhs[4]);
    
    /* Create a matrix for the return argument */
    plhs[0] = mxCreateDoubleMatrix(ndim,mdim,mxREAL);
    dCum = mxGetPr(plhs[0]);
        
    /* Do the actual computations in a subroutine */
   mkCumProd(d,dCum,pwt,mdim,ndim,mlim); 
    
}


