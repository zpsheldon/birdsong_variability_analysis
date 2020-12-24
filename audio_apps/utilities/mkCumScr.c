/*=================================================================
 *
 * mkCumScr.c
 *	        cumulative score for dynamic time warping
 *          paths.
 *
 * The calling syntax is:
 *
 *		D = mkCumScr(d,pwt,mdim,ndim,mlim)
 *      d: raw distance matrix
 *      pwt: path weightings
 *      mdim/ndim: d rows/cols
 *      mlim,nlim: limits defining Sakoe-Chiba band (Rabiner & Juang 1993)
 *
 *      values outside the band are 0 by default
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2000 The MathWorks, Inc.
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"

void mkCumScr(
        double	*d,
        double   *dCum,
        double  *pwt,
        double mdim,
        double ndim,
        double mlim,
        double nlim
        ) {
    
    double	val1, val2, val3, dimfct;
    int n, m, nlim2, lm1, lm2, p11, p12, p23, p21, p22, p32;
    
    
    dimfct = mdim / ndim;
    
    /* Populate first row/column with values from distance matrix */
    for (m = 1; m <= mlim; m++) {
        p11 = (m-1)*ndim;
        *(dCum+p11) = *(d+p11);
    }
    
    nlim2 = ceil(1/dimfct)+nlim < ndim? ceil(1/dimfct)+nlim:ndim;
    
    for (n = 1; n <= nlim2; n++) {
        p11 = n-1;
        *(dCum+p11) = *(d+p11);
    }
    
    
    /* point (2,2) in cumulative is just d(1,1)+d(2,2)*/
    p11 = 1+ndim;
    *(dCum+p11) = *(d+p11) + *(d);
    
    
    /* 2nd row and column need to be calculated separately because
     * some paths from main algorithm point outside matrix*/
    n = 2;
    lm1 = n-mlim;
    lm2 = mlim+n-1;
    lm1 = lm1>3? lm1:3;
    lm2 = lm2<mdim? lm2:mdim;
    for (m = lm1; m <= lm2; m++) {
        
        p11 = n+(m-1)*ndim-1;
        
        p12 = n+(m-2)*ndim-1;
        p23 = n-1+(m-3)*ndim-1;
        
        p22 = n-1+(m-2)*ndim-1;
        
        val1 = *(pwt) * *(d+p11) + *(dCum+p22);
        val3 = *(pwt+2) * (.5 * *(d+p11) + .25 * *(d+p12) + .25 * *(d+p22)) + *(dCum+p23);
        
        *(dCum+p11) = val1>val3? val1 : val3;
        
    }
    
    m = 2;
    nlim2 = ceil(2/dimfct)+nlim < ndim? ceil(2/dimfct)+nlim:ndim;
    for (m = lm1; m <= lm2; m++) {
        
        p11 = n+(m-1)*ndim-1;
        
        p21 = n-1+(m-1)*ndim-1;
        p32 = n-2+(m-2)*ndim-1;
        
        p22 = n-1+(m-2)*ndim-1;
        
        val1 = *(pwt) * *(d+p11) + *(dCum+p22);
        val2 = *(pwt+1) * (.5 * *(d+p11) + .25 * *(d+p21) + .25 * *(d+p22)) + *(dCum+p32);
        
        *(dCum+p11) = val1>val2? val1 : val2;
        
    }
    
    
    /* main algorithm*/
    for (n = 3; n <= ndim; n++) {
        
        lm1 = ceil((n-1)*dimfct)-mlim;
        lm2 = ceil((n-1)*dimfct)+mlim;
        
        lm1 = lm1>3? lm1:3;
        lm2 = lm2<mdim? lm2:mdim;
        
        for (m = lm1; m <= lm2; m++) {
            
            p11 = n+(m-1)*ndim-1;
            
            p21 = n-1+(m-1)*ndim-1;
            p32 = n-2+(m-2)*ndim-1;
            
            p12 = n+(m-2)*ndim-1;
            p23 = n-1+(m-3)*ndim-1;
            
            p22 = n-1+(m-2)*ndim-1;
            
            val1 = *(pwt) * *(d+p11) + *(dCum+p22);
            val2 = *(pwt+1) * (.5 * *(d+p11) + .25 * *(d+p21) + .25 * *(d+p22)) + *(dCum+p32);
            val3 = *(pwt+2) * (.5 * *(d+p11) + .25 * *(d+p12) + .25 * *(d+p22)) + *(dCum+p23);
            
            *(dCum+p11) = val1>val2? (val1>val3? val1:val3) : (val2>val3? val2:val3);
            
        }
        
    }
    
}

void mexFunction( int nlhs, mxArray*plhs[],
        int nrhs, const mxArray*prhs[])
        
{
    double *dCum, *d;
    double *pwt;
    int mdim, ndim, mlim, nlim;
    
    d = mxGetPr(prhs[0]);
    pwt = mxGetPr(prhs[1]);
    
    mdim = mxGetScalar(prhs[2]);
    ndim = mxGetScalar(prhs[3]);
    mlim = mxGetScalar(prhs[4]);
    nlim = mxGetScalar(prhs[5]);
    
    /* Create a matrix for the return argument */
    plhs[0] = mxCreateDoubleMatrix(ndim, mdim, mxREAL);
    dCum = mxGetPr(plhs[0]);
    
    /* Do the actual computations in a subroutine */
    mkCumScr(d, dCum, pwt, mdim, ndim, mlim, nlim);
    
}


