/*=================================================================
 *
 * BAUM-WELCH ALGORITHM FOR ESTIMATION OF POSTERIORS, ETC IN HMM
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void baumwelch(
        double   *Apr,
        double  *Bpr,
        double *ppr,
        int K,
        int N,
        int D,
        double *alphapr,
        double *betapr,
        double *cpr,
        double *gammacovpr
        ) {
    
    /* E is DXK, X is DXN, A is KXK, p is KX1, alpha and beta KXN, gammacov KXK */
    int n, k, k2, d;
    
    
    /* INITIALIZE ALPHA */
    for (k = 0; k < K; k++) {
        *(alphapr + k) = *(Bpr + k) * *(ppr + k);
        *cpr += *(alphapr + k);
    }
    
    *cpr = 1 / *cpr;
    
    for (k = 0; k < K; k++) {
        *(alphapr + k) *= *cpr;
    }
    
    
    /* ALPHA RECURSION */
    
    for (n = 1; n < N; n++) {
        for (k = 0; k < K; k++) {
            
            for (k2 = 0; k2 < K; k2++) {
                *(alphapr + k + n*K) += *(alphapr + k2 + (n-1)*K) * *(Apr + k2 + k*K);
            }
            
            *(alphapr + k + n*K) *= *(Bpr + n*D + k);
            *(cpr + n) += *(alphapr + k + n*K);
        }
        
        *(cpr + n) = 1 / *(cpr + n);
        
        for (k = 0; k < K; k++) {
            *(alphapr + k + n*K) *= *(cpr + n);
        }
        
    }
    
    
    /* BETA INITIALIZATION */
    
    for (k = 0; k < K; k++) {
        *(betapr + k + (N-1)*K) = *(cpr + N - 1);
    }
    
    
    /* BETA RECURSION */
    
    for (n = N-1; n > 0; n--) {
        for (k = 0; k < K; k++) {
            
            for (k2 = 0; k2 < K; k2++) {
                *(betapr + k + (n-1)*K) += *(Bpr + n*K + k2) * *(betapr + k2 + n*K) * *(Apr + k + k2*K);
            }
            
            *(betapr + k + (n-1)*K) *= *(cpr + n - 1);
        }
        
        
    }
    
    
    /* GAMMA 2ND MOMENT*/
    
    for (n = 0; n < N - 1; n++) {
        for (k = 0; k < K; k++) {
            for (k2 = 0; k2 < K; k2++) {
                *(gammacovpr + k + k2*K) += *(Bpr + (n+1)*K + k2) * *(alphapr + k + n*K) * *(Apr + k + k2*K) * *(betapr + k2 + (n+1)*K);
            }
            
        }
        
        
    }
    
    
}

void mexFunction( int nlhs, mxArray*plhs[],
        int nrhs, const mxArray*prhs[])
        
{
    double *Apr, *Bpr, *ppr, *alphapr, *betapr, *gammacovpr, *cpr;
    int N, D, K;
    
    Apr = mxGetPr(prhs[0]);
    Bpr = mxGetPr(prhs[1]);
    ppr = mxGetPr(prhs[2]);
    
    
    K = mxGetN(prhs[0]);
    N = mxGetN(prhs[1]);
    D = mxGetM(prhs[1]);
    
    /*argnm = mxGetN(prhs);*/
    
    
    /* Create matrices for the return argument */
    plhs[0] = mxCreateDoubleMatrix(K, N, mxREAL);
    alphapr = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(K, N, mxREAL);
    betapr = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1, N, mxREAL);
    cpr = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(K, K, mxREAL);
    gammacovpr = mxGetPr(plhs[3]);
    
    
    /* Do the actual computations in a subroutine */
    baumwelch(Apr, Bpr, ppr, K, N, D, alphapr, betapr, cpr, gammacovpr);
    
}


