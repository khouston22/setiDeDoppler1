
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"
//#include "matrix.h"
#include "kiss_fft.h"

#if !defined(PI)
#define PI (3.14159265358979323846)
#endif

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

// function [x_out,fs_out,f_bin,t_out] = run_pfb_Lfx(x_in,coef,fs_in,n_out,pLf,qLf);
  /* 
  /*  Polyphase Filter Bank with critically sampled output and
  /*  rational (integer ratio Lf = pLf/qLf) bin overlap
  /*
  /* Function performs coefficient summation
  /* FFTs done separately 
  /* 
  /* inputs
  /* 
  /* xr_in,xi_in    n_in          input waveform (complex float)
  /* coef           M*N_tap       coefficients (real float)
  /* xcoefr,xcoefi  M*N_tap*pLf   coefficients work area (complex float)
  /* workr,worki    n_freq x 1    pfb work vector (complex float)
  /* n_in           scalar        number of PFB output samples per subfilter
  /* n_out          scalar        number of PFB output samples per subfilter
  /* pLf,qLf        scalar        Lf = pLf/qLf filter bank overlap factor
  /*                              =1 for standard critically sampled FB ("1x")
  /*                              =2 for 50% overlapped frequency bins ("2x")
  /*                              rational number: Lf = pLf/qLf
  /*                                                                                                                                                    
  /* outputs
  /* 
  /* xr_out,xi_out  n_freq*n_out    output waveform (complex)
  /*
  /* All arrays are treated as vectors
  /* Complex vectors are treated as two real vectors (r,i)
  /*
  /* pfb_Lf1 performs polyphase filter summation
  /* pfb_Lf2 performs FFTs and data shuffling
  */
/*  ======================================================================  */

void pfb_Lf1(float *xr_sum, float *xi_sum, long n_out, 
             float *xr_in, float *xi_in, long n_in, 
             float *coef, float *xcoefr, float *xcoefi, 
             float *workr, float *worki, 
             long M, long N_tap, long pLf, long qLf)
{
long i1,i_tap,k,ofs1,ofs2;
long n_freq = M*pLf/qLf;
long N_fft = M/qLf;
double arg;
float c1r,c1i,c2r,c2i;
double Lf = (double) pLf/(double) qLf;

/* coefficient matrix including phase offsets */
/* effectively 3D matrix:  M x N_tap x pLf */
/* xcoef(:,:,i1+1) = (coef/M).*(exp(-1j*2*pi*i1*[0:M-1]'/n_freq)*exp(-1j*2*pi*[0:N_tap-1]*i1/Lf)); */

for (i1=0; i1<pLf; i1++) {
  ofs2 = i1*M*N_tap;
  ofs1 = 0;
  if (i1==0) {
    for (i_tap=0; i_tap<N_tap; i_tap++) {
      for (k=0; k<M; k++) {
        xcoefr[ofs2]   = coef[ofs1++]/M;
        xcoefi[ofs2++] = 0.;
      }
    }
  } else {
    for (i_tap=0; i_tap<N_tap; i_tap++) {
      arg = -2.0*PI*i_tap*i1/Lf;
      c2r = (float) cos(arg);
      c2i = (float) sin(arg);
      for (k=0; k<M; k++) {
        arg = -2.0*PI*k*i1/n_freq;
        c1r = coef[ofs1]  /M * (float) cos(arg);
        c1i = coef[ofs1++]/M * (float) sin(arg);
        xcoefr[ofs2]   = c1r*c2r - c1i*c2i;
        xcoefi[ofs2++] = c1r*c2i + c1i*c2r;
      }
    }
  }
}

/* do the filter bank summations */

  float xr,xi,xcr,xci;
  long i_out,ofs_out,ofs_in,ofs_c,ofs_w;
  long n_in_max = M*(n_in/M)-M;
          
  for (i_out=0; i_out<n_out; i_out++) {
    for (i1=0; i1<pLf; i1++) {
      memset(workr,0,M*sizeof(float));
      memset(worki,0,M*sizeof(float));
      ofs_out = i_out*n_freq + i1*N_fft;
      ofs_in  = i_out*M;
      ofs_c = i1*M*N_tap;
      for (i_tap=0; i_tap<N_tap; i_tap++) {
        if (ofs_in>(n_in_max)) break;  // be sure not to go over end of x_in
        for (k=0; k<M; k++) {
          xr = xr_in[ofs_in];
          xi = xi_in[ofs_in++];
          xcr = xcoefr[ofs_c];
          xci = xcoefi[ofs_c++];
          workr[k] += xr*xcr - xi*xci;
          worki[k] += xr*xci + xi*xcr; 
        }
      }
      memcpy(&xr_sum[ofs_out],&workr[0],N_fft*sizeof(float));
      memcpy(&xi_sum[ofs_out],&worki[0],N_fft*sizeof(float));
      for (long j1=1; j1<qLf; j1++) {
        ofs_w = j1*N_fft;
        for (k=0; k<N_fft; k++) {
          xr_sum[ofs_out+k] += workr[ofs_w];
          xi_sum[ofs_out+k] += worki[ofs_w++];
        }
      }
    }
  }
  
  return;
}
/* PFB FFT and output data shuffle function 
/* function is in-place, so xr_sum xi_sum input and xr_out xi_out are same array
/* workr worki needs to be n_freq long and can be same array as for pfb_Lf1
/* This version use kiss_fft for FFTs
*/
void pfb_Lf2(float *xr_out, float *xi_out, float *workr, float *worki, 
             kiss_fft_cpx *fft_in, kiss_fft_cpx *fft_out,
             kiss_fft_cfg st,
             long n_out, long M, long pLf, long qLf)
{
  long ofs_out,i_out,i1,k;
  long n_freq = M*pLf/qLf;
  long N_fft = M/qLf;
  
  for (i_out=0; i_out<n_out; i_out++) {
    /* copy summation data for one time instant (one column) */
    ofs_out = i_out*n_freq;
    memcpy(&workr[0],&xr_out[ofs_out],n_freq*sizeof(float));
    memcpy(&worki[0],&xi_out[ofs_out],n_freq*sizeof(float));

    /* copy input to output vector before running in-place FFT */
    for (i1=0; i1<pLf; i1++) {
      for (k=0; k<N_fft; k++) {
        fft_in[k].r = workr[i1*N_fft+k];
        fft_in[k].i = worki[i1*N_fft+k];
      }

      kiss_fft(st,fft_in,fft_out);

      /* copy FFT output to work, and perform FFT shift */

      for (k=0; k<N_fft/2; k++) {
        workr[i1*N_fft+k] = fft_out[k+N_fft/2].r;
        worki[i1*N_fft+k] = fft_out[k+N_fft/2].i;
        workr[i1*N_fft+k+N_fft/2] = fft_out[k].r;
        worki[i1*N_fft+k+N_fft/2] = fft_out[k].i;
      }
    }
    
    /* copy work output back to final output, while interleaving */

    for (i1=0; i1<pLf; i1++) {
      for (k=0; k<N_fft; k++) {
        xr_out[ofs_out+k*pLf+i1] = workr[i1*N_fft+k];
        xi_out[ofs_out+k*pLf+i1] = worki[i1*N_fft+k];
      }
    }
  }
  
  return;
}
