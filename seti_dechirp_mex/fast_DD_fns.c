
#include <math.h>
#include <string.h>
#include <stdlib.h>
//#include "mex.h"
//#include "matrix.h"

#if !defined(MAX)
#define	MAX(A, B)	(((A) > (B)) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	(((A) < (B)) ? (A) : (B))
#endif

void fastDD_stage1(float *det_DD1, float *xx_ext, long m_min0, long m_max0, 
               long Nt, long N0, long n_freq1, long n_freq2, 
               long LfLr_ratio, long alg_ID, long dd0);
void fastDD_stage_i(float *det_DD_out, float *det_DD_in, long i_stage,
                   long m_min0, long m_max0,long Nt, long N0, long n_freq2, 
                   long alg_ID, long dd0);
void m_limit(long *m_min_i,long *m_max_i,long m_min0,long m_max0,long Ni,long N0);


/*  ======================================================================  */
/*  These functions perform a De-Doppler sum on a spectrogram. It assumes  */
/*  the arrangement of data stream is, all points in first spectra, all     */
/*  points in second spectra, etc...  Data are summed across time           */
/*                     Original version: K. Houston 2023                    */
/*  ======================================================================  */

void fastDD_stage1(float *det_DD1, float *xx_ext, long m_min0, long m_max0, 
               long Nt, long N0, long n_freq1, long n_freq2, 
               long LfLr_ratio, long alg_ID, long dd0)
{
  /* 
  /*  fastDD "Slow" first stage
  /* 
  /* det_DD1[(Nr0*n_group)*(n_freq2+2*dd0)]  output array
  /* xx_ext[Nt*(n_freq1+2*dd0)]              input spectrogram (zero padded)
  /* 
  /* Note that det_DD1 array is assumed zeroed with calloc() call
  */
  long ig,it0,m0,it1,it,m,mm,mm0,kk,kk0,t_g_ofs,m_g_ofs;
  long n_freq1_ext = n_freq1 + 2*dd0;
  long n_freq2_ext = n_freq2 + 2*dd0;
  long n_group = Nt/N0;
  long Nr0 = m_max0 - m_min0 + 1;
  long *int_bin_offset;
  long xx_ofs,dd_ofs;

  /*  ======================================================================  */

  /* calculate bin offsets */
    
//   if (alg_ID>=1) {
    int_bin_offset = (long *) malloc((size_t)Nr0*N0*sizeof(long));
    for (mm=0; mm<Nr0; mm++) {
      m = m_min0 + mm;
      for (it=0; it<N0; it++) {
        int_bin_offset[mm*N0+it] = (long) round((double)LfLr_ratio*m*it/(N0));
      }
    }
//   } else {  // Taylor 1st stage
//     int_bin_offset = (long *) malloc((size_t)Nr0*N0*sizeof(long));
//     for (mm=0; mm<Nr0; mm++) {
//       m = m_min0 + mm;
//       for (it=0; it<N0; it++) {
//         int_bin_offset[mm*N0+it] = (long) round((double)LfLr_ratio*m*it/(N0-1));
//       }
//     }
//   }

  for (ig=0; ig<n_group; ig++) {
    t_g_ofs = ig*N0;
    m_g_ofs = ig*Nr0;
    for (mm=0; mm<Nr0; mm++) {
      mm0 = m_g_ofs+mm;
      for (it=0; it<N0; it++) {
        it1 = t_g_ofs+it;
        // input SG offset equivalent to xx_ext[mm0][dd0+int_bin_offset[mm][it]]
        xx_ofs = it1*n_freq1_ext + dd0 + int_bin_offset[mm*N0+it];
        // output DD offset equivalent to det_DD1[mm0][dd0]
        dd_ofs = mm0*n_freq2_ext + dd0; 
        
        if (it==0) {
          for (kk=0; kk<n_freq2; kk++) {
            det_DD1[dd_ofs++] = xx_ext[xx_ofs];
            xx_ofs += LfLr_ratio;
          }
        } else {
          for (kk=0; kk<n_freq2; kk++) {
            det_DD1[dd_ofs++] += xx_ext[xx_ofs];
            xx_ofs += LfLr_ratio;
          }
        }
      }
    }
  }
  free(int_bin_offset);
  return;
}

void fastDD_stage_i(float *det_DD_out, float *det_DD_in, long i_stage,
                   long m_min0, long m_max0,long Nt, long N0, long n_freq2, 
                   long alg_ID, long dd0)
{
  /* 
  /*  "Fast" DeDoppler stage
  /* 
  /* det_DD_out[(Nr1*n_group1<=Nr)*(n_freq2+2*dd0)]  output array - stage i_stage
  /* det_DD_in[(Nr0*n_group0<=Nr)*(n_freq2+2*dd0)]  input array - stage i_stage-1
  /* 1 <= alg_ID <= 6
  */
  long ig,m1,kk,i1,i2;
  long k_ofs0,k_ofs1,k_ofs2;
  long m_g_ofs10,m_g_ofs11,m_g_ofs2;
  long dd_ofs10,dd_ofs11,dd_ofs2;
  long dd_ofs10a,dd_ofs11a,dd_ofs2a;
  long dd_ofs10b,dd_ofs11b,dd_ofs2b;
  long dd_ofs10c,dd_ofs11c,dd_ofs2c;
  long n_freq2_ext = n_freq2 + 2*dd0;
  long N_out = N0<<(i_stage-1);  //   N2 = N0* 2.^(i_stage-1);
  long N_in = N_out/2;  
  long n_group = Nt/N_out;
  long m_min_out,m_max_out,m_min_in,m_max_in;

  /*  ======================================================================  */

  m_max_out = (long) round(m_max0*N_out/N0);
  m_min_out = (long) round(m_min0*N_out/N0);
  m_max_in = (long) round(m_max0*N_in/N0);
  m_min_in = (long) round(m_min0*N_in/N0);

  long Nr_out = m_max_out - m_min_out + 1;  // number of drift rates in current stage
  long Nr_in  = m_max_in - m_min_in + 1;    // number of drift rates in previous stage

  for (ig=0; ig<n_group; ig++) {
    m_g_ofs10 = ig*2*Nr_in;
    m_g_ofs11 = m_g_ofs10 + Nr_in;
    m_g_ofs2  = ig*Nr_out;

    for (i1=0; i1<Nr_in; i1++) {
      m1 = i1 + m_min_in;
      k_ofs0 = dd0 + m1;
      k_ofs1 = dd0 + m1 + 1;
      k_ofs2 = dd0 + m1 + 2;

      // "even" output points
      i2 = 2*i1;
      // input DD offset equivalent to 
      // det_DD_in[m_g_ofs10+i1][dd0] & det_DD_in[m_g_ofs11+i1][k_ofs0]
      dd_ofs10 = (m_g_ofs10+i1)*n_freq2_ext + dd0; 
      dd_ofs11 = (m_g_ofs11+i1)*n_freq2_ext + k_ofs0; 
      
      // output DD offset equivalent to det_DD_out[m_g_ofs2+i2][dd0]
      dd_ofs2 = (m_g_ofs2+i2)*n_freq2_ext + dd0; 

      for (kk=0; kk<n_freq2; kk++) {
        det_DD_out[dd_ofs2++] = det_DD_in[dd_ofs10++] + det_DD_in[dd_ofs11++];
      }
      
      i2 = 2*i1+1;
      if (i2<Nr_out) {
        // "odd" output points
        // output DD offset equivalent to det_DD_out[m_g_ofs2+i2][dd0]
        dd_ofs2 = (m_g_ofs2+i2)*n_freq2_ext + dd0; 
        if (alg_ID==1) {
          // point-by-point max of two dog-leg sums
          dd_ofs10a = (m_g_ofs10+i1  )*n_freq2_ext + dd0;
          dd_ofs10b = (m_g_ofs10+i1+1)*n_freq2_ext + dd0; 
          dd_ofs11a = (m_g_ofs11+i1+1)*n_freq2_ext + k_ofs0; 
          dd_ofs11b = (m_g_ofs11+i1  )*n_freq2_ext + k_ofs1; 

          for (kk=0; kk<n_freq2; kk++) {
            det_DD_out[dd_ofs2++] = 
                    fmax((det_DD_in[dd_ofs10a++] + det_DD_in[dd_ofs11a++]),
                         (det_DD_in[dd_ofs10b++] + det_DD_in[dd_ofs11b++]));
          }
        } else if (alg_ID==2) {
          // cubic interpolation
          dd_ofs10a = (m_g_ofs10+i1  )*n_freq2_ext + dd0; 
          dd_ofs10b = (m_g_ofs10+i1+1)*n_freq2_ext + dd0; 
          dd_ofs10c = (m_g_ofs10+i1+2)*n_freq2_ext + dd0; 
          dd_ofs11a = (m_g_ofs11+i1+1)*n_freq2_ext + k_ofs0; 
          dd_ofs11b = (m_g_ofs11+i1  )*n_freq2_ext + k_ofs1; 
          dd_ofs11c = (m_g_ofs11+i1-1)*n_freq2_ext + k_ofs2; 
          if (i1==0) {
            for (kk=0; kk<n_freq2; kk++) {
              det_DD_out[dd_ofs2++] = .375*det_DD_in[dd_ofs10a++] + 
                                      .750*det_DD_in[dd_ofs10b++] + 
                                     -.125*det_DD_in[dd_ofs10c++] + 
                                      .375*det_DD_in[dd_ofs11a++] +
                                      .750*det_DD_in[dd_ofs11b++];
            }
          } else if (i1==Nr_in-2) {
            for (kk=0; kk<n_freq2; kk++) {
              det_DD_out[dd_ofs2++] = .375*det_DD_in[dd_ofs10a++] + 
                                      .750*det_DD_in[dd_ofs10b++] + 
                                      .375*det_DD_in[dd_ofs11a++] +
                                      .750*det_DD_in[dd_ofs11b++] +
                                     -.125*det_DD_in[dd_ofs11c++];
            }
          } else {
            for (kk=0; kk<n_freq2; kk++) {
              det_DD_out[dd_ofs2++] = .375*det_DD_in[dd_ofs10a++] + 
                                      .750*det_DD_in[dd_ofs10b++] + 
                                     -.125*det_DD_in[dd_ofs10c++] + 
                                      .375*det_DD_in[dd_ofs11a++] +
                                      .750*det_DD_in[dd_ofs11b++] +
                                     -.125*det_DD_in[dd_ofs11c++];
            }
          }
        } else if (alg_ID>=3) {
          dd_ofs2 = (m_g_ofs2+i2)*n_freq2_ext + dd0; 
          if (alg_ID==3) {
            // lower dog leg
            dd_ofs10 = (m_g_ofs10+i1  )*n_freq2_ext + dd0; 
            dd_ofs11 = (m_g_ofs11+i1+1)*n_freq2_ext + k_ofs0; 
          } else if (alg_ID==4) {
            // upper dog leg
            dd_ofs10 = (m_g_ofs10+i1+1)*n_freq2_ext + dd0; 
            dd_ofs11 = (m_g_ofs11+i1  )*n_freq2_ext + k_ofs1; 
          } else if (alg_ID==5) {
            // upper leg then lower
            dd_ofs10 = (m_g_ofs10+i1+1)*n_freq2_ext + dd0; 
            dd_ofs11 = (m_g_ofs11+i1+1)*n_freq2_ext + k_ofs0; 
         } else if (alg_ID==6) {
            // lower leg then upper
            dd_ofs10 = (m_g_ofs10+i1+1)*n_freq2_ext + dd0; 
            dd_ofs11 = (m_g_ofs11+i1  )*n_freq2_ext + k_ofs1; 
          }
          for (kk=0; kk<n_freq2; kk++) {
            det_DD_out[dd_ofs2++] = det_DD_in[dd_ofs10++] + det_DD_in[dd_ofs11++];
          }
        }
      }
    }
  }

  return;
}

