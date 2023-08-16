/* ======================================================================  */
/* This is a mex function to run a polyphase filter bank (PFB) on an input */
/* signal, reference run_pfb_Lfx.m                                         */
// function [x_out,fs_out,f_bin,t_out] = run_pfb_Lfx(x_in,coef,fs_in,n_out,pLf,qLf);
// %
// % function to run Polyphase Filter Bank with input x_in
// % Generates interpolated output bins
// %
// % inputs
// %
// % x_in         n_in x 1   input waveform
// % coef         MxN_tap    coefficients
// % fs_in        1x1        PFB input sample rate, Hz
// % n_out        1x1        number of PFB output samples per subfilter
// % pLf,qLf      1x1        Lf = pLf/qLf filter bank overlap factor
// %                         =1 for standard critically sampled FB ("1x")
// %                         =2 for 50% overlapped frequency bins ("2x")
// %                         rational number: Lf = pLf/qLf                                                                                                                                                    
// % outputs
// %
// % x_out        n_freq x n_out  output waveform
// % fs_out       1x1             filter bank output sample rate, Hz = fs_in/M
// % f_bin        n_freq x 1      bin center freq values
// % t_out        n_out x 1       output time adjusted for PFB delay N_tap/2/fs_out
// %


#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"
#include "kiss_fft.h"

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#if !defined(DEBUG)
#define DEBUG 0
#endif

#if !defined(PI)
#define PI (3.14159265358979323846)
#endif

void pfb_Lf1(float *xr_out, float *xi_out, long n_out, 
             float *xr_in, float *xi_in, long n_in, 
             float *coef, float *xcoefr, float *xcoefi, 
             float *workr, float *worki, 
             long M, long N_tap, long pLf, long qLf);

void pfb_Lf2(float *xr_out, float *xi_out, float *workr, float *worki, 
             kiss_fft_cpx *fft_in, kiss_fft_cpx *fft_out,
             kiss_fft_cfg st,
             long n_out, long M, long pLf, long qLf);

// function [x_out,fs_out,f_bin,t_out] = pfb_Lfx(x_in,coef,fs_in,n_out,pLf,qLf);

/* Input Arguments to mex call */

#define	X_IN	 prhs[0]
#define	COEF	 prhs[1]
#define	FS_IN  prhs[2]
#define	N_OUT  prhs[3]
#define	PLF    prhs[4]
#define	QLF    prhs[5]
#define	DO_FFT prhs[6]

/* Output Arguments */

#define	X_OUT	 plhs[0]
#define	FS_OUT plhs[1]
#define	F_BIN	 plhs[2]
#define T_OUT	 plhs[3]
#define XCOEF	 plhs[4]

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  /* Check for proper number of arguments */

  if (nrhs != 7) { 
    mexErrMsgIdAndTxt( "MATLAB:pfb_Lfx_mex:invalidNumInputs",
              "Six input arguments required."); 
  } else if (nlhs > 5) {
    mexErrMsgIdAndTxt( "MATLAB:pfb_Lfx_mex:maxlhs",
              "Too many output arguments."); 
  } 

  size_t x_in_rows, x_in_cols;
  x_in_rows = mxGetM(X_IN); 
  x_in_cols = mxGetN(X_IN); 
  long n_in = MAX(x_in_rows,x_in_cols);

  if (!mxIsSingle(X_IN)) { 
    mexErrMsgIdAndTxt( "MATLAB:pfb_Lfx_mex:invalidX",
              "pfb_Lfx_mex requires that X_IN be a vector of type SINGLE."); 
  } 
    
  double fs_in = mxGetScalar(FS_IN);
  long   n_out  = (long) round(mxGetScalar(N_OUT));
  long   pLf    = (long) round(mxGetScalar(PLF));
  long   qLf    = (long) round(mxGetScalar(QLF));
  long   do_fft = (long) round(mxGetScalar(DO_FFT));
  long   M      = mxGetM(COEF); 
  long   N_tap  = mxGetN(COEF); 
  long   n_in_min  = (n_out + N_tap-1)*M; 
  
  double  Lf = pLf/qLf;
  long n_freq = M*pLf/qLf;
  long N_fft = M/qLf; 
  double fs_out = fs_in/M;

  #if DEBUG
    mexPrintf("\nChecking inputs:\n");
    mexPrintf("x_in: %ld x %ld, n_in=%ld, n_out=%ld\n",x_in_rows,x_in_cols,n_in,n_out);
    mexPrintf("Lf=%.2f,pLf=%ld,qLf=%ld,M=%ld,N_tap=%ld\n",Lf,pLf,qLf,M,N_tap);
    mexPrintf("fs_in=%.2f,fs_out=%.2f\n\n",fs_in,fs_out);
  #endif
    
  /* Get pointer to x_in data */ 

  float  *xr, *xi, *coef, *xcoefr, * xcoefi, *workr, * worki, *zr, *zi;
  xr = (float *) mxGetData(X_IN);
  xi = (float *) mxGetImagData(X_IN);
  coef = (float *) mxGetData(COEF);
  
  #if DEBUG
    mexPrintf("Checking x_in values\n");
    mexPrintf("x_in[0]=%.1f %.1f,x_in[1]=%.1f %.1f,x_in[2]=%.1f %.1f,x_in[3]=%.1f %.1f\n",
       xr[0],xi[0],xr[1],xi[1],xr[2],xi[2],xr[3],xi[3]);
  #endif
    
  if (M*pLf!=n_freq*qLf) { 
    mexErrMsgIdAndTxt( "MATLAB:pfb_Lfx_mex:invalid_M_qLf",
              "pfb_Lfx_mex requires that M*pLf/qLf must be integer"); 
  } 

  /* Allocate output data */ 
    
  X_OUT = mxCreateNumericMatrix( n_freq,n_out,mxSINGLE_CLASS,mxCOMPLEX);
  zr = (float *) mxGetData(X_OUT);
  zi = (float *) mxGetImagData(X_OUT);
  
  /* Allocate work data */ 

  xcoefr = (float *) mxMalloc((mwSize) M*N_tap*pLf*sizeof(float));
  xcoefi = (float *) mxMalloc((mwSize) M*N_tap*pLf*sizeof(float));
  workr = (float *) mxMalloc((mwSize) n_freq*sizeof(float));
  worki = (float *) mxMalloc((mwSize) n_freq*sizeof(float));

  /* Call pfb summartion function */ 

  pfb_Lf1(zr,zi,n_out,xr,xi,n_in,coef,xcoefr,xcoefi,workr,worki,M,N_tap,pLf,qLf);
  
  /* Do FFTs for pfb */ 

  if (do_fft) {    
    int k;
    int nfft[32];
    int ndims = 1;
    int isinverse=0;
    kiss_fft_cpx *fft_in;
    kiss_fft_cpx *fft_out;
    int real = 0;
    nfft[0] = N_fft;
    int nbytes = sizeof(kiss_fft_cpx);
    for (k=0;k<ndims;++k)
      nbytes *= nfft[k];
    #if DEBUG
      mexPrintf("sizeof(kiss_fft_cpx)=%ld, nbytes=%ld\n",sizeof(kiss_fft_cpx),nbytes);
    #endif

    fft_in =(kiss_fft_cpx*)KISS_FFT_MALLOC(nbytes);
    fft_out=(kiss_fft_cpx*)KISS_FFT_MALLOC(nbytes);
    kiss_fft_cfg st = kiss_fft_alloc( nfft[0] ,isinverse ,0,0);
    
    pfb_Lf2(zr,zi,workr,worki,fft_in,fft_out,st,n_out,M,pLf,qLf);
    
    free(st);
    free(fft_in); free(fft_out);
    kiss_fft_cleanup();
  }

  /* Output additional arrays */ 

  if (nlhs > 1) {
    FS_OUT = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS, mxREAL);
    double *fs_out_p;
    fs_out_p = (double *) mxGetPr(FS_OUT); 
    *fs_out_p = fs_out;
  }

  double *f_bin;
  if (nlhs > 2) {
    F_BIN = mxCreateNumericMatrix(n_freq,1,mxDOUBLE_CLASS, mxREAL);
    f_bin = (double *) mxGetPr(F_BIN); 
    for (long k=0; k<n_freq; k++) {
      f_bin[k] = (k-(n_freq/2))*fs_in/n_freq;
    }
  }

  double *t_out;
  if (nlhs > 3) {
    T_OUT = mxCreateNumericMatrix(1,n_out,mxDOUBLE_CLASS, mxREAL);
    t_out = (double *) mxGetPr(T_OUT); 
    for (long n=0; n<n_out; n++) {
      t_out[n] = (n+N_tap/2)/fs_out;
    }
  }

  mwSize n_dims,dims[3];
  n_dims = 3;
  dims[0] = M;
  dims[1] = N_tap;
  dims[2] = pLf;
  float *xcoefr_p,*xcoefi_p;

  if (nlhs > 4) {
    XCOEF = mxCreateNumericArray( n_dims,dims,mxSINGLE_CLASS,mxCOMPLEX);
    xcoefr_p = (float *) mxGetData(XCOEF);
    xcoefi_p = (float *) mxGetImagData(XCOEF);
    memcpy(xcoefr_p,xcoefr, (size_t) M*N_tap*pLf*sizeof(float));
    memcpy(xcoefi_p,xcoefi, (size_t) M*N_tap*pLf*sizeof(float));
  }

  mxFree(xcoefr);
  mxFree(xcoefi);
  mxFree(workr);
  mxFree(worki);
  
  return;
}

