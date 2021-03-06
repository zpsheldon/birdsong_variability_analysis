/*=================================================================

    2-D CROSS-CORRELATION THAT IS ONLY MATCHES ALONG COLUMNS. ARGUMENTS ARE: MATRIX 1, MATRIX 2, MAXIMUM LAG TIME, AND WHETHER
    THE MATCH IS NORMALIZED BY THE NUMBER OF POINTS CALCULATED IN THE MATCH AT A POSITION.
 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void xcorr21(
		   double	*y,
           double   *x,
           double  *c,
           double *lags,
           int ycols,
           int xcols,
           int ccols,
           int rows,
           int maxlag,
           double *nvc
		   )
{

    int xind1,yind1,clen,col,row,n,counter;
    double p;
    
    for (n = 1; n <= ccols; n++) {
        
        *(lags+n-1) = n-maxlag-1;
        
        xind1 = *(lags+n-1)*rows;
        xind1 = xind1 > 0 ? xind1:0;
        
        yind1 = (maxlag-n+1)*rows;
        yind1 = yind1 > 0 ? yind1:0;
        
        clen = *(lags+n-1) < 0 ? ycols + *(lags+n-1):ycols;
        clen = xcols - *(lags+n-1) > clen ? clen: xcols - *(lags+n-1);

        p = 0;
        counter = 0;
        
        for (col = 1; col <= clen; col++) {
            for (row = 1; row <= rows; row++) {

              /*p += pow(*(y+yind1+counter) - *(x+xind1+counter),2)*/;
              p += (*(y+yind1+counter) - *(x+xind1+counter)) * (*(y+yind1+counter) - *(x+xind1+counter));
              counter++;
        
            }
            
        }
        
        p = p / counter;
        
        *(c+n-1) = p;
        *(nvc+n-1) = counter;

    }
        
}

void mexFunction( int nlhs, mxArray*plhs[], 
		  int nrhs, const mxArray*prhs[])
     
{ 
    double *x,*y,*c,*maxlag,*lags,*nvc;
    int xcols,ycols,ccols,yrows,xrows;
    
    y = mxGetPr(prhs[0]);
    x = mxGetPr(prhs[1]);
    maxlag = mxGetPr(prhs[2]);
    
    ycols = mxGetN(prhs[0]);
    xcols = mxGetN(prhs[1]);
    yrows = mxGetM(prhs[0]);
    xrows = mxGetM(prhs[1]);
    
    /*argnm = mxGetN(prhs);*/
    
    if (yrows != xrows) {
    mexErrMsgTxt("x and y must have same row number.");
  }
    
    ccols = xcols-ycols+1+2*maxlag[0];
    
    /* Create a matrix for the return argument */
    plhs[0] = mxCreateDoubleMatrix(1,ccols,mxREAL);
    c = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,ccols,mxREAL);
    lags = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1,ccols,mxREAL);
    nvc = mxGetPr(plhs[2]);
    
        
    /* Do the actual computations in a subroutine */
   xcorr21(y,x,c,lags,ycols,xcols,ccols,xrows,maxlag[0],nvc); 
    
}


