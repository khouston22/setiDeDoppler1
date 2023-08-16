/* ======================================================================  */
/* This is a mex function to Taylor-tree-sum a data stream, 
/* using an alternative algorithm fastDD, reference fastDD3.m */

// function [det_DD,freq2,df_dt_list,m_list] = ...
//                  fastDD3(xx,freq1,T_line,Lf,Lr,Nt,N0,df_dt_min,df_dt_max,alg_ID);
// %
// % Function to run energy detection with input xx and 
// % a range of Doppler drift rates
// % 
// % Outputs
// %
// % inputs
// %
// % xx           n_freq1 x n_time  mag squared time waveforms for each frequency bin (power)
// % freq1        n_freq1 x 1       PFB bin center freq values
// % T_line       1 x 1             time sec for each line in spectrogram
// % Lf           1 x 1             overlap factor in spectrogram (power of 2)
// % Lr           1 x 1             upsample for drift rate increment
// % Nt            1 x 1             last stage number of time samples (Nt/N0 is power of 2)
// % N0           1 x 1             first stage number of time samples
// % df_dt_min    1 x 1             minimum drift rate to be computed >= -fs^2 (approx)
// % df_dt_max    1 x 1             maximum drift rate to be computed <= fs^2 (approx)
// % alg_ID       1 x 1             ID for algorithm variant, =1 baseline
// %                                =0 extended Taylor DD
// %
// % outputs
// %
// % det_DD       n_freq2 x n_m     integrated energy array
// % freq2        n_freq2 x 1       det_DD bin center freq values
// % df_dt_list   n_m x 1           det_DD drift rate Hz/sec array
// % m_list       n_m x 1           drift rate index array
// %
// % where
// %
// % n_freq2 = n_freq1/max(Lf/Lr,1)
// % n_m = depends on df_dt_min, df_dt_max
// %
// % Notes
// %
// % class of det_DD will be same as class of input xx, e.g. 'double' or 'single'
// %
// % Normally n_time will equal Nt.  If n_time<Nt, input will be zero padded and
// % only n_time lines will be summed. If n_time>Nt input will be truncated and
// % only Nt integrations will be performed.
// %
// % Define bin_bw = bin bandwidth in each pfb bin = fs = sampling rate at
// % bin output
// %
// % Input in xx will have frequency oversampling factor Lf, so 
// %    delta_freq1 = bin_bw/Lf = freq1(2) - freq1(1)
// % Output in det_DD does critical frequency sampling factor, so 
// %    delta_freq2 = bin_bw/Lr = freq2(2) - freq2(1)
// % Note Lf >= Lr, and Lf must be an integer multiple of Lr
// %
// % Output will also have a df_dt increment
// %    delta_df_dt = df_dt_list(2) - df_dt_list(1),
// % where delta_df_dt = fs/(Lr*n_time*T_line)
// %
// % Normally T_line will equal 1/fs.  Sometimes the spectrogram lines will be
// % averaged N_preDD times prior to input, so that T_line = N_preDD/fs.
// % However, this is generally not recommended (though accommodated here),
// % because the frequency drifts over fs^2/N_preDD will suffer significant
// % attenuation in the averaging process.
// % 

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#if !defined(DEBUG)
#define DEBUG 0
#endif

void taylor_stage1(float *det_DD1, float *xx_ext, long m_min0, long m_max0, 
               long Nt, long N0, long n_freq1, long n_freq2, 
               long LfLr_ratio, long dd0);
void taylor_stage_i(float *det_DD_out, float *det_DD_in, long i_stage,
                   long m_min0, long m_max0,long Nt, long N0, long n_freq2, 
                   long dd0);
void m_limit(long *m_min_i,long *m_max_i,long m_min0,long m_max0,long Ni,long N0);
void fastDD_stage1(float *det_DD1, float *xx_ext, long m_min0, long m_max0, 
               long Nt, long N0, long n_freq1, long n_freq2, 
               long LfLr_ratio, long alg_ID, long dd0);
void fastDD_stage_i(float *det_DD_out, float *det_DD_in, long i_stage,
                   long m_min0, long m_max0,long Nt, long N0, long n_freq2, 
                   long alg_ID, long dd0);

// function [det_DD,freq2,df_dt_list,m_list] = ...
//                  fastDD3(xx,freq1,T_line,Lf,Lr,Nt,N0,df_dt_min,df_dt_max,alg_ID);

/* Input Arguments to mex call */

#define	XX_IN	prhs[0]
#define	FREQ_IN	prhs[1]
#define	TLine_IN	prhs[2]
#define	Lf_IN	prhs[3]
#define	Lr_IN	prhs[4]
#define	Nt_IN	prhs[5]
#define	N0_IN	prhs[6]
#define	dFdTmin_IN	prhs[7]
#define	dFdTmax_IN	prhs[8]
#define	alg_ID_IN	prhs[9]

/* Output Arguments */

#define	DetDD_OUT	plhs[0]
#define	FREQ_OUT	plhs[1]
#define	dFdT_OUT	plhs[2]
#define	mIndex_OUT	plhs[3]

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  /* Check for proper number of arguments */

  if (nrhs != 10) { 
    mexErrMsgIdAndTxt( "MATLAB:fast_DD:invalidNumInputs",
              "Eight input arguments required."); 
  } else if (nlhs > 4) {
    mexErrMsgIdAndTxt( "MATLAB:fast_DD:maxlhs",
              "Too many output arguments."); 
  } 

  long n_freq1 = mxGetM(XX_IN); 
  long Nt_ = mxGetN(XX_IN);

  if (!mxIsSingle(XX_IN) || (n_freq1==1) || (Nt_==1)) { 
    mexErrMsgIdAndTxt( "MATLAB:fast_DD:invalidXX",
              "fast_DD requires that XX_IN be a matrix of type SINGLE."); 
  } 

  long n_freq1_r = mxGetM(FREQ_IN); 
  long n_freq1_c = mxGetN(FREQ_IN);
  long n_freq1_ = MAX(n_freq1_r,n_freq1_c);
  if (!mxIsDouble(FREQ_IN) || (n_freq1!=n_freq1_)) { 
    mexPrintf("n_freq1=%ld vs %ld\n",n_freq1,n_freq1_);
    mexErrMsgIdAndTxt( "MATLAB:fast_DD:invalidFreq1",
              "fast_DD requires that freq1 be a vector of type Double."); 
  }  
    
  double T_line = mxGetScalar(TLine_IN);
  double Lf = mxGetScalar(Lf_IN);
  double Lr = mxGetScalar(Lr_IN);
  long   Nt = (long) round(mxGetScalar(Nt_IN));
  long   N0 = (long) round(mxGetScalar(N0_IN));
  double df_dt_min = mxGetScalar(dFdTmin_IN);
  double df_dt_max = mxGetScalar(dFdTmax_IN);
  long   alg_ID = (long) round(mxGetScalar(alg_ID_IN));

  long LfLr_ratio = (long) round(Lf/Lr);
  long n_freq2 = n_freq1/LfLr_ratio;

  #if DEBUG
    mexPrintf("Checking inputs,n_freq1=%ld,n_freq2=%ld,Nt=%ld,N0=%ld,Lf=%.2f,Lr=%.2f,LfLr_ratio=%ld,alg_ID=%ld\n",
         n_freq1,n_freq2,Nt,N0,Lf,Lr,LfLr_ratio,alg_ID);
  #endif

  if (fabs(LfLr_ratio*Lr-Lf)>1e-5) { 
    mexErrMsgIdAndTxt( "MATLAB:fast_DD:invalidLfLr_ratio",
              "fast_DD requires that Lf/Lr be an integer"); 
  }  
   
  if ((alg_ID>6) || (alg_ID<0)) { 
    mexErrMsgIdAndTxt( "MATLAB:fast_DD:invalidalg_ID",
              "fast_DD requires 0 <= alg_ID <=6"); 
  }  
   
  /* Copy frequency data into input matrix */ 

  mxArray *freq1_p;
  double *freq1;
  
  freq1_p = mxCreateNumericArray(mxGetNumberOfDimensions(FREQ_IN),
            mxGetDimensions(FREQ_IN),mxDOUBLE_CLASS, mxREAL);
  freq1 = (double *) mxGetPr(freq1_p); 
  memcpy(freq1, mxGetPr(FREQ_IN), (size_t) n_freq1*mxGetElementSize(FREQ_IN));
  
  /*
  /* compute constants
  */

  double df_bin = freq1[1]-freq1[0];  // may be less than 1/Ts if bins are overlapped
  double fs = df_bin*Lf;
  long n_stage = (long) round(log2((double)(Nt/N0)))+1;
  
  #if DEBUG
    mexPrintf("Checking freq parameters\n");
    mexPrintf("freq1[0]=%.2f,freq1[1]=%.2f,df_bin=%.2f,Lf=%.2f,fs=%.2f,n_stage=%ld\n",
       freq1[0],freq1[1],df_bin,Lf,fs,n_stage);
  #endif
  
  long m_min0,m_max0,m_min,m_max;
  
  if (alg_ID>=1) {
    m_max0 = (long) ceil((double)df_dt_max*T_line/fs*(Lr*N0));
    m_min0 = (long) floor((double)df_dt_min*T_line/fs*(Lr*N0));
    m_max = (long) round(m_max0*Nt/N0);
    m_min = (long) round(m_min0*Nt/N0);
  } else {
    m_max0 = (long) ceil((double)df_dt_max*T_line/fs*(Lr*N0-1));
    m_min0 = (long) floor((double)df_dt_min*T_line/fs*(Lr*N0-1));
    m_limit(&m_min,&m_max,m_min0,m_max0,Nt,N0);
  }
  long Nr = m_max - m_min + 1; 
  long Nr_ext = MAX((m_max0 - m_min0 + 1)*Nt/N0,Nr); 

  long dd0 = MAX(abs(m_max),abs(m_min));
  dd0 = (long) pow(2,ceil(log2((double)dd0))+1.);
  
  #if DEBUG
    mexPrintf("Checking drift rate limits\n");
    mexPrintf("df_dt_min=%.2f,df_dt_max=%.2f,T_line=%.2f,fs=%.2f,Lr=%.2f,N0=%ld\n",
       df_dt_min,df_dt_max,T_line,fs,Lr,N0);
    mexPrintf("m_min0=%ld,m_max0=%ld,m_min=%ld,m_max=%ld,Nr=%ld,Nr_ext=%ld,n_stage=%ld,dd0=%ld\n",
       m_min0,m_max0,m_min,m_max,Nr,Nr_ext,n_stage,dd0);
  #endif

  /* Copy spectrogram data into input matrix, with zero padding */ 

  float *xx,*xx_ext;
  mwSize n_dims,dims[2];
  dims[0] = Nt;
  dims[1] = n_freq1;

  xx = (float *) mxGetData(XX_IN);

  #if DEBUG
    mexPrintf("Checking SG values\n");
    mexPrintf("xx[0]=%.2f,xx[1]=%.2f,xx[2]=%.2f,xx[3]=%.2f\n",
       xx[0],xx[1],xx[2],xx[3]);
  #endif
  
  /* Allocate and zero-pad spectrogram data */ 

  long n_freq1_ext = n_freq1 + 2*dd0;
    
  xx_ext = (float *) mxMalloc((mwSize) Nt*n_freq1_ext*sizeof(float)); // no init to zero

  for (long it=0; it<Nt_; it++) {
    memset(&xx_ext[it*n_freq1_ext],0,dd0*sizeof(float));
    memcpy(&xx_ext[it*n_freq1_ext + dd0], &xx[n_freq1*it], (size_t) n_freq1*sizeof(float));
    memset(&xx_ext[it*n_freq1_ext+dd0+n_freq1],0,dd0*sizeof(float));
  }
  
  #if DEBUG
    mexPrintf("Checking zero pad SG values\n");
    mexPrintf("zp %.2f %.2f, xx_ext %.2f %.2f %.2f %.2f\n",
       xx_ext[0],xx_ext[1], xx_ext[dd0+0],xx_ext[dd0+1],xx_ext[dd0+2],xx_ext[dd0+3]);
  #endif
  
  /* allocate work array for output */
  
  float *det_DD_work[2];
  long n_freq2_ext = n_freq2 + 2*dd0;
  long i_in=1;
  long i_out=0;
  long i_stage=1;

  det_DD_work[0] = (float *) mxCalloc(Nr_ext*n_freq2_ext,sizeof(float));
  det_DD_work[1] = (float *) mxCalloc(Nr_ext*n_freq2_ext,sizeof(float));
  
  if (alg_ID>=1) {
    /* first stage processing - ffastDD algorithm */

    fastDD_stage1(det_DD_work[0],xx_ext,m_min0,m_max0,Nt,N0,n_freq1,n_freq2,LfLr_ratio,alg_ID,dd0);

    /* i-th stage processing */

    for (i_stage=2; i_stage<=n_stage; i_stage++) {
      i_in = i_stage % 2;
      i_out = (i_stage+1) % 2;
      #if DEBUG
        mexPrintf("Begin stage %ld,i_in=%ld,i_out=%ld\n",i_stage,i_in,i_out);
      #endif
      fastDD_stage_i(det_DD_work[i_out],det_DD_work[i_in],i_stage,m_min0,m_max0,
                    Nt,N0,n_freq2,alg_ID,dd0);
    }
  } else {
    /* first stage processing - Taylor alg_ID=0 */

    taylor_stage1(det_DD_work[0],xx_ext,m_min0,m_max0,Nt,N0,n_freq1,n_freq2,LfLr_ratio,dd0);

    /* i-th stage processing */

    for (i_stage=2; i_stage<=n_stage; i_stage++) {
      i_in = i_stage % 2;
      i_out = (i_stage+1) % 2;
      #if DEBUG
        mexPrintf("Begin stage %ld,i_in=%ld,i_out=%ld\n",i_stage,i_in,i_out);
      #endif
      taylor_stage_i(det_DD_work[i_out],det_DD_work[i_in],i_stage,m_min0,m_max0,
                    Nt,N0,n_freq2,dd0);
    }
  }

  /* Create output matrices for the return arguments */ 

  float *detDD_out,*freq_out,*df_dt_list,*m_list;
  DetDD_OUT = mxCreateNumericMatrix(n_freq2,Nr,mxSINGLE_CLASS, mxREAL);
  detDD_out = (float *) mxGetData(DetDD_OUT);
 
  for (long mm=0; mm<Nr; mm++) {
    for (long i_freq=0; i_freq<n_freq2; i_freq++) {
      detDD_out[mm*n_freq2+i_freq] = det_DD_work[i_out][mm*n_freq2_ext+dd0+i_freq]/Nt;
    }
  }
  
  if (nlhs > 1) {
    if (n_freq1_r>n_freq1_c) {
      FREQ_OUT = mxCreateNumericMatrix(n_freq2,1,mxSINGLE_CLASS,mxREAL);
    } else {
      FREQ_OUT = mxCreateNumericMatrix(1,n_freq2,mxSINGLE_CLASS,mxREAL);
    }
    freq_out = (float *) mxGetData(FREQ_OUT); 
    for (long if2=0; if2<n_freq2; if2++) {
      freq_out[if2] = freq1[LfLr_ratio*if2];
    }
  }

  if (nlhs > 2) {
    dFdT_OUT = mxCreateNumericMatrix(1,Nr,mxSINGLE_CLASS,mxREAL);
    df_dt_list = (float *) mxGetData(dFdT_OUT); 
    for (long m=m_min; m<=m_max; m++) {
      long mm = m - m_min;
      df_dt_list[mm] = m*fs/T_line/(Nt*Lr);
    }
  }

  if (nlhs > 3) {
    mIndex_OUT = mxCreateNumericMatrix(1,Nr,mxSINGLE_CLASS,mxREAL);
    m_list = (float *) mxGetData(mIndex_OUT);
    for (long m=m_min; m<=m_max; m++) {
      long mm = m - m_min;
      m_list[mm] = m;
    }
  }

  mxFree(xx_ext);
  mxFree(det_DD_work[0]);
  mxFree(det_DD_work[1]);
  
  return;
}
