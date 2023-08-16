/*  ======================================================================  */
/*  This is a function to Taylor-tree-sum a data stream. Adapted from    */
/* https://github.com/UCBerkeleySETI/gbt_seti/blob/master/src/rawdopplersearch.c */

#include <math.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"


void taylor_flt(float outbuf[], long int mlen, long int nchn);
long int bitrev(long int inval,long int nbits);

/* Input Arguments */

#define	INBUF	prhs[0]
#define	MLEN_IN	prhs[1]
#define	NCHN_IN	prhs[2]
#define	OUTBUF	plhs[0]


/* Output Arguments */

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif


/*  ======================================================================  */
/*  This function bit-reverses the given value "inval" with the number of   */
/*  bits, "nbits".    ----  R. Ramachandran, 10-Nov-97, nfra.               */
/*  ======================================================================  */

long int bitrev(long int inval,long int nbits)
{
     long int     ifact,k,i,ibitr;

     if(nbits <= 1)
     {
          ibitr = inval;
     }
     else
     {
          ifact = 1;
          for (i=1; i<(nbits); ++i)
               ifact  *= 2;
          k     = inval;
          ibitr = (1 & k) * ifact;

          for (i=2; i < (nbits+1); i++)
          {
               k     /= 2;
               ifact /= 2;
               ibitr += (1 & k) * ifact;
          }
     }
     return ibitr;
}

/*  ======================================================================  */
/*  This is a function to Taylor-tree-sum a data stream. It assumes that    */
/*  the arrangement of data stream is, all points in first spectra, all     */
/*  points in second spectra, etc...  Data are summed across time           */
/*                     Original version: R. Ramachandran, 07-Nov-97, nfra.  */
/*                     Modified 2011 A. Siemion float/64 bit addressing     */
/*  outbuf[]       : input array (float), replaced by dedispersed data  */
/*                   at the output                                          */
/*  mlen           : dimension of outbuf[] (long int)                            */
/*  nchn           : number of frequency channels (long int)                     */
/*                                                                          */
/*  ======================================================================  */

void taylor_flt(float outbuf[], long int mlen, long int nchn)
{
  float itemp;
  long int   nsamp,npts,ndat1,nstages,nmem,nmem2,nsec1,nfin, i;
  long int   istages,isec,ipair,ioff1,i1,i2,koff,ndelay,ndelay2;
  long int   bitrev(long int, long int);

  /*  ======================================================================  */

  nsamp   = ((mlen/nchn) - (2*nchn));
  npts    = (nsamp + nchn);
  ndat1   = (nsamp + 2 * nchn);
  
  //nstages = (int)(log((float)nchn) / 0.6931471 + 0.5);
  nstages = (long int) log2((double)nchn);
  nmem    = 1;


  for (istages=0; istages<nstages; istages++) {
    nmem  *= 2;
    nsec1  = (nchn/nmem);
    nmem2  = (nmem - 2);

    for (isec=0; isec<nsec1; isec++) {
      ndelay = -1;
      koff   = (isec * nmem);

      for (ipair=0; ipair<(nmem2+1); ipair += 2) {
        

        ioff1   = (bitrev(ipair,istages+1)+koff) * ndat1;
        i2      = (bitrev(ipair+1,istages+1) + koff) * ndat1;
        ndelay++;
        ndelay2 = (ndelay + 1);
        nfin    = (npts + ioff1);
        for (i1=ioff1; i1<nfin; i1++) {
          itemp      = (outbuf[i1] + outbuf[i2+ndelay]);
          outbuf[i2] = (outbuf[i1] + outbuf[i2+ndelay2]);
          outbuf[i1] = itemp;
          i2++;

        }
      }
    }
  }

  return;
}


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    size_t m,n; 
  
    mxArray *work_p;
    float *work; 
    float *out; 
    double dmlen,dnchn;
    long mlen,nchn,n_freq,n_time,nbits,ibrev; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 3) { 
	    mexErrMsgIdAndTxt( "MATLAB:taylor_flt:invalidNumInputs",
                "Three input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "MATLAB:taylor_flt:maxlhs",
                "Too many output arguments."); 
    } 
    
    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */ 
    
    m = mxGetM(INBUF); 
    n = mxGetN(INBUF);

    if (!mxIsSingle(INBUF) || (MAX(m,n) == 1) || (MIN(m,n) != 1) 
          || (MIN(m,n) != 1)) { 
	    mexErrMsgIdAndTxt( "MATLAB:taylor_flt:invalidINBUF",
                "taylor_flt requires that INBUF be a vector of type SINGLE."); 
    } 
    if (!mxIsDouble(MLEN_IN) || !mxIsScalar(MLEN_IN)) { 
	    mexErrMsgIdAndTxt( "MATLAB:taylor_flt:invalidMLEN",
                "taylor_flt requires that MLEN be a scalar of type Double."); 
    } 
    if (!mxIsDouble(NCHN_IN) || !mxIsScalar(NCHN_IN)) { 
	    mexErrMsgIdAndTxt( "MATLAB:taylor_flt:invalidNCHN",
                "taylor_flt requires that NCHN be a scalar of type Double."); 
    } 
    
    /* Create a matrix for the return argument */ 
    
    work_p = mxCreateNumericArray(mxGetNumberOfDimensions(INBUF), 
              mxGetDimensions(INBUF),mxSINGLE_CLASS, mxREAL);
    OUTBUF = mxCreateNumericArray(mxGetNumberOfDimensions(INBUF), 
              mxGetDimensions(INBUF),mxSINGLE_CLASS, mxREAL);

    /* Assign pointers to the various parameters */ 

    work = (float *) mxGetPr(work_p); 
    memcpy(work, mxGetPr(INBUF), (size_t) m*n*mxGetElementSize(INBUF));

    memcpy(&dmlen, mxGetPr(MLEN_IN), mxGetElementSize(MLEN_IN));
    memcpy(&dnchn, mxGetPr(NCHN_IN), mxGetElementSize(NCHN_IN));

    mlen = (long) dmlen;
    nchn = (long) dnchn;
    n_freq = mlen/nchn;
    n_time = nchn;
    
    if (0) mexPrintf("Before call,mlen=%ld,nchn=%ld,n_time=%ld,n_freq=%ld\n",
            mlen,nchn,n_time,n_freq);
    if (0) {
      for (long i=0;i<MIN(1024,mlen);i++) {
        mexPrintf("%3.0f ",work[i]);
        if (i%16 == 15) {
          mexPrintf("\n");
        }
      }
      mexPrintf("\n");
    }
    
//     for (long nbits=1; nbits<8; nbits++) {
//       mexPrintf("bitrev test,bits=%ld\n",nbits);
//       for (i=0;i<pow(2.0,nbits);i++) {
//         mexPrintf("%3ld ",bitrev(i,nbits));
//         if (i%16 == 15) {
//           mexPrintf("\n");
//         }
//       }
//       mexPrintf("\n");
//     }
       
    /* Do the computations */
    
    taylor_flt(work, mlen, nchn);

    if (0) {
      mexPrintf("After call,mlen=%ld,nchn=%ld\n",mlen,nchn);
      for (long i=0;i<MIN(1024,mlen);i++) {
        mexPrintf("%3.0f ",work[i]);
        if (i%16 == 15) {
          mexPrintf("\n");
        }
      }
      mexPrintf("\n");
    }

    /* bit-reverse the time output */
    
    out = (float *) mxGetPr(OUTBUF); 
    nbits = (long) log2((double)n_time);

    for (long i=0; i<n_time; i++) {
      ibrev = bitrev(i, nbits); 
      //mexPrintf("Bit Reversal,nbits=%ld,i=%ld,ibrev=%ld\n",nbits,i,ibrev);
      memcpy(&out[ibrev*n_freq], &work[i*n_freq], (size_t) n_freq*sizeof(float));
    }

    if (0) {
      mexPrintf("After call,mlen=%ld,nchn=%ld\n",mlen,nchn);
      for (long i=0;i<MIN(1024,mlen);i++) {
        mexPrintf("%3.0f ",out[i]);
        if (i%16 == 15) {
          mexPrintf("\n");
        }
      }
      mexPrintf("\n");
    }

    return;
    
}
