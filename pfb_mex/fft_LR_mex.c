/* ======================================================================  */
/* This is a mex function to perform an FFT                                */
// function [x_out] = fft_LR_mex(x_in,N);

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"
#include "fft_LR.h"

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#if !defined(DEBUG)
#define DEBUG 0
#endif

/* Input Arguments to mex call */

#define	X_IN	prhs[0]
#define	N_IN	prhs[1]

/* Output Arguments */

#define	X_OUT	plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  /* Check for proper number of arguments */

  if (nrhs != 2) { 
    mexErrMsgIdAndTxt( "MATLAB:fft_LR_mex:invalidNumInputs",
              "Eight input arguments required."); 
  } else if (nlhs > 1) {
    mexErrMsgIdAndTxt( "MATLAB:fft_LR_mex:maxlhs",
              "Too many output arguments."); 
  } 

  size_t rows, cols;
  rows = mxGetM(X_IN); 
  cols = mxGetN(X_IN); 

  double dN;
  long N;

  memcpy(&dN, mxGetPr(N_IN), (size_t) mxGetElementSize(N_IN));

  N = (long) round(dN);

  if (!mxIsSingle(X_IN)) { 
    mexErrMsgIdAndTxt( "MATLAB:fft_LR_mex:invalidXX",
              "fft_LR_mex requires that X_IN be a matrix of type SINGLE."); 
  } 
    
  #if DEBUG
    mexPrintf("Checking inputs,x_in: %ld x %ld vs N=%ld\n\n",rows,cols,N);
  #endif

  if (rows!=N) { 
    mexPrintf("x_in: %ld x %ld vs N=%ld\n",rows,cols,N);
    mexErrMsgIdAndTxt( "MATLAB:fft_LR_mex:invalidN",
              "fft_LR_mex requires that X_IN be a N x n_col matrix"); 
  } 
       
  /* Get pointer to x_in data */ 

  float  *xr, *xi, *zr, *zi;
  xr = (float *) mxGetData(X_IN);
  xi = (float *) mxGetImagData(X_IN);
  
  #if DEBUG
    mexPrintf("Checking x_in values\n");
    mexPrintf("x_in[0]=%.1f %.1f,x_in[1]=%.1f %.1f,x_in[2]=%.1f %.1f,x_in[3]=%.1f %.1f\n",
       xr[0],xi[0],xr[1],xi[1],xr[2],xi[2],xr[3],xi[3]);
  #endif
  
  /* Allocate output data */ 
  mwSize n_dims,dims[2];
  n_dims = 2;
  dims[0] = rows;
  dims[1] = cols;

  X_OUT = mxCreateNumericArray( n_dims,dims,mxSINGLE_CLASS,mxCOMPLEX);
  zr = (float *) mxGetData(X_OUT);
  zi = (float *) mxGetImagData(X_OUT);
  
  /* repeat for each column */
  
  for (long i_col=0; i_col<cols; i_col++) {
    /* copy input to output vector before running in-place FFT */

    memcpy(&zr[i_col*N],&xr[i_col*N], (size_t) N*sizeof(float));
    memcpy(&zi[i_col*N],&xi[i_col*N], (size_t) N*sizeof(float));

    /* run the FFT function - data is in-place */

    fft(&zr[i_col*N], &zi[i_col*N], N);
  }

  return;
}

