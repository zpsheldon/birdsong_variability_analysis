/*=================================================================
 *
 * MATCHTEMPLATE.C	
 * SIMLAR TO XCORR21, EXCEPT THAT LAG TIME IS RESTRICTED TO ONSET OF SPECTROGRAM TO BE MATCHED
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2000 The MathWorks, Inc.
 *
 *=================================================================*/
/* $Revision: 1.10 $ */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void xseqc(
        double	*y,
        double   *x,
        double  *c,
        int ycols,
        int xcols,
        int ccols
        ) {
    
    int xcol, ycol, counter;
    double p;
    
    for (xcol = 0; xcol < ccols; xcol++) {
        
        p = 0;
        counter = 0;
        
        for (ycol = 0; ycol < ycols; ycol++) {
            
            p += fabs(*(y+ycol) - *(x+xcol+ycol));
            counter++;
            
        }
        
        if (p==0){*(c+xcol) = 1;}
        
    }
    
    
    
}

void mexFunction( int nlhs, mxArray*plhs[],
		  int nrhs, const mxArray*prhs[])
     
{ 
    double *x,*y,*c;
    int xcols,ycols,ccols;
    
    y = mxGetPr(prhs[0]);
    x = mxGetPr(prhs[1]);
    
    ycols = mxGetN(prhs[0]);
    xcols = mxGetN(prhs[1]);
    
    ccols = xcols-ycols+1;
    
    /* Create a matrix for the return argument */
    plhs[0] = mxCreateDoubleMatrix(1,ccols,mxREAL);
    c = mxGetPr(plhs[0]);
        
    /* Do the actual computations in a subroutine */
   xseqc(y,x,c,ycols,xcols,ccols); 
    
}


