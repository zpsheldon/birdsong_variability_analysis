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

void matchTemplate(
		   double	*y,
           double   *x,
           double  *c,
           int ycols,
           int xcols,
           int ccols,
           int rows
		   )
{

    int xind1,col,row,n,counter;
    double p;
    
    for (n = 1; n <= ccols; n++) {
        
        xind1 = (n-1)*rows;
        
        p = 0;
        counter = 0;
        
        for (col = 1; col <= ycols; col++) {
            for (row = 1; row <= rows; row++) {

              /*p += fabs(*(y+counter) - *(x+xind1+counter));*/
              p += pow((*(y+counter) - *(x+xind1+counter)),2);
                
              counter++;
        
            }
            
        }
        
        *(c+n-1) = counter/p;
         
    }
        
}

void mexFunction( int nlhs, mxArray*plhs[], 
		  int nrhs, const mxArray*prhs[])
     
{ 
    double *x,*y,*c;
    int xcols,ycols,ccols,yrows,xrows;
    
    y = mxGetPr(prhs[0]);
    x = mxGetPr(prhs[1]);
    
    ycols = mxGetN(prhs[0]);
    xcols = mxGetN(prhs[1]);
    yrows = mxGetM(prhs[0]);
    xrows = mxGetM(prhs[1]);
    
    if (yrows != xrows) {
    mexErrMsgTxt("x and y must have same row number.");
  }
    
    ccols = xcols-ycols+1;
    
    /* Create a matrix for the return argument */
    plhs[0] = mxCreateDoubleMatrix(1,ccols,mxREAL);
    c = mxGetPr(plhs[0]);
        
    /* Do the actual computations in a subroutine */
   matchTemplate(y,x,c,ycols,xcols,ccols,xrows); 
    
}


