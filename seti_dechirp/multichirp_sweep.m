%
% DeDoppler of multiplexed chirps in one input waveform
% with one taylorDD1 or fastDD call
%

%
% generate run PFB on input signal
%

fprintf('PFB%.0fx %3.0f  ',Lf,n_chirp);
tstart = tic;
[s_out,freq1,t_out,coef] = gen_chirp_pfb_matrix2(s_in,fs,Nt_pfb,M,N_tap,window_name,Lf);        
mean_t_cpu_pfb(i_N0,i_case) = toc(tstart);
fprintf(' %.1f sec\n',mean_t_cpu_pfb(i_N0,i_case));

t_max = n_in/fs_in;
t_in = [0:n_in-1]/fs_in;

%
% loop over snr
%

for i_snr=1:n_snr
  
  snr_db = snr_db_list(i_snr);

  fprintf('Case %.0f of %.0f, %.0f Chirps, SNR %.0f of %.0f SNR=%.1f dB\n',...
          i_case,n_case,n_chirp,i_snr,n_snr,snr_db);

  if isinf(snr_db)
    x_out = A_list(i_snr)*s_out;
  else
    x_out = A_list(i_snr)*s_out + sigma_list(i_snr)*noise_out;
  end

  xx = abs(x_out).^2;

  %
  % do pre-DD averaging
  %

  if (N_pre_DD>1)
    [xx,t_out] = pre_DD_avg(xx,t_out,N_pre_DD);
  end
  
  %
  % do RFI suppression
  %

  if rfi_proc_enable
    xx_rfi = xx;
    [xx,bin_power_in,bin_power_out] = rfi_agc1(xx,N_rfi_avg_f,N_rfi_avg_t,rfi_weight_exp,...
                      rfi_bin_limit,rfi_bin_replace,bin_noise_magsq0);
  end
  
  %xx_db = 10*log10(max(1e-4,xx));

  %
  % compute pfb signal envelopes
  %

  xx_ssg_cat_db = [];

  for i_chirp = 1:n_chirp

    df_extract = f_incr_offset/2;
    [xx_ssg,f_ssg] = shift_sg(xx,freq,t_out,f1_all(i_chirp),df_dt_all(i_chirp),df_extract);

    xx_ssg_cat_db = [xx_ssg_cat_db ; 10*log10(xx_ssg)];

    taylor_best_mean = mean(max(xx_ssg));

    taylor_best_db(i_chirp,i_snr,i_case) = 10*log10(taylor_best_mean);
    bin_noise_db(i_chirp,i_snr,i_case) = bin_noise_db0;
    
    if isinf(snr_db)
      rel_snr_best_db(i_chirp,i_case) = taylor_best_db(i_chirp,i_snr,i_case) - bin_noise_db0;
    end
  end

  %
  % "fast" energy detector = "TaylorDD"
  % 

  fprintf('Taylor-%.0f v%.0f   ',N0_taylor,ver_DD);

  tstart = tic;

  if (ver_DD==4)
    [det_taylor,freq_taylor,df_dt_list_taylor,m_list_t] = ...
                 taylor_DD_mex(xx,freq,T_line,Lf,Lr,N_DD,N0_taylor,df_dt_min,df_dt_max);
  elseif (ver_DD==3)
    [det_taylor,freq_taylor,df_dt_list_taylor,m_list_t] = ...
                 taylorDD3(xx,freq,T_line,Lf,Lr,N_DD,N0_taylor,df_dt_min,df_dt_max);
  end

  t_cpu_taylor(1:n_chirp,i_snr,i_case) = toc(tstart);
  fprintf('%.1f sec, ',t_cpu_taylor(1,i_snr,i_case));

  if (i_snr==1)&&(i_case==1)
    det_taylor_max_vs_f_db = [];
    det_taylor_SE_vs_f_db = [];
    det_taylor_max_vs_df_dt_db = [];
    det_taylor_SE_vs_df_dt_db = [];
    freq_taylor1 = freq_taylor;
    df_dt_list_taylor1 = df_dt_list_taylor;
  end

  det_taylor_max_vs_f_db{i_snr,i_case} = 10*log10(max(.01,max(det_taylor,[],2)));
  det_taylor_SE_vs_f_db{i_snr,i_case} = ...
                    10*log10((max(.01,max(det_taylor,[],2)-bin_noise_magsq))/mdl_taylor);

  det_taylor_max_vs_df_dt_db{i_snr,i_case} = 10*log10(max(.01,max(det_taylor,[],1)));
  det_taylor_SE_vs_df_dt_db{i_snr,i_case} = ...
                    10*log10((max(.01,max(det_taylor,[],1)-bin_noise_magsq))/mdl_taylor);

  %
  % extract peaks from taylor detection array
  %

  dfreq_taylor = freq_taylor(2)-freq_taylor(1);
  nbf2 = round(f_incr_offset/dfreq_taylor/2);
  nbf1 = 2*nbf2+1;
  nmb2 = 10;
  nmb1 = 2*nmb2+1;

%   det_taylor_x_all = NaN(nbf1,nmb1,n_chirp);
%   det_taylor_max_vs_f = NaN(nbf1,n_chirp);
%   det_taylor_max_vs_m = NaN(nmb1,n_chirp);

  for i_chirp = 1:n_chirp

    if (df_dt_all(i_chirp)>max(df_dt_list_taylor))
      continue;
    end

    % postage stamp of det array
    [temp,i_freq_taylor] = min(abs(freq_taylor-f1_all(i_chirp)));
    [temp,i_m0] = min(abs(df_dt_list_taylor-df_dt_all(i_chirp)));
    ii_m = max(1,i_m0-nmb2):min(size(det_taylor,2),i_m0+nmb2);

    det_taylor_x = det_taylor(i_freq_taylor+[-nbf2:nbf2],ii_m);

    freq_taylor_x = freq_taylor(i_freq_taylor+[-nbf2:nbf2]);
    df_dt_taylor_x = df_dt_list_taylor(ii_m);
    m_list_taylor_x = m_list_t(ii_m);

    [det_max1,i1,i2] = max2d(det_taylor_x);
    taylor_db(i_chirp,i_snr,i_case) = 10*log10(det_max1);
    taylor_f_est(i_chirp,i_snr,i_case) = freq_taylor_x(i1);
    taylor_m_est(i_chirp,i_snr,i_case) = m_list_taylor_x(i2);
    taylor_df_dt_est(i_chirp,i_snr,i_case) = df_dt_taylor_x(i2);

    if isinf(snr_db)
      rel_snr_taylor_db(i_chirp,i_case) = taylor_db(i_chirp,i_snr,i_case) - bin_noise_db0;
    else
      SE_taylor_db(i_chirp,i_snr,i_case) = 10*log10((det_max1 - bin_noise_magsq)/mdl_taylor);
    end
    
    taylor_f_error(i_chirp,i_snr,i_case) = freq_taylor_x(i1)-f1_all(i_chirp);
    taylor_df_dt_error(i_chirp,i_snr,i_case) = df_dt_taylor_x(i2)-df_dt_all(i_chirp);
  end

  %
  % "fast" De-Doppler = "fastDD"
  % 

  fprintf('fastDD-%.0f v%.0f   ',N0,ver_DD);

  tstart = tic;

  if (ver_DD==4)
    [det_fastDD,freq_fastDD,df_dt_list_fastDD,m_list_fastDD] = ...
                   fast_DD_mex(xx,freq,T_line,Lf,Lr,N_DD,N0,df_dt_min,df_dt_max,fastDD_alg_ID);
  elseif (ver_DD==3)
    [det_fastDD,freq_fastDD,df_dt_list_fastDD,m_list_fastDD] = ...
                   fastDD3(xx,freq,T_line,Lf,Lr,N_DD,N0,df_dt_min,df_dt_max,fastDD_alg_ID);
  end

  t_cpu_fastDD(1:n_chirp,i_snr,i_case) = toc(tstart);
  fprintf('%.1f sec, alg ID=%.0f\n',t_cpu_fastDD(1,i_snr,i_case),fastDD_alg_ID);
 
  if (i_snr==1)&&(i_case==1)
    det_fastDD_max_vs_f_db = [];
    det_fastDD_SE_vs_f_db = [];
    det_fastDD_max_vs_df_dt_db = [];
    det_fastDD_SE_vs_df_dt_db = [];
  end

  det_fastDD_max_vs_f_db{i_snr,i_case} = 10*log10(max(.01,max(det_fastDD,[],2)));
  det_fastDD_SE_vs_f_db{i_snr,i_case} = ...
                    10*log10((max(.01,max(det_fastDD,[],2)-bin_noise_magsq))/mdl_fastDD);
  det_fastDD_max_vs_df_dt_db{i_snr,i_case} = 10*log10(max(.01,max(det_fastDD,[],1)));
  det_fastDD_SE_vs_df_dt_db{i_snr,i_case} = ...
                    10*log10((max(.01,max(det_fastDD,[],1)-bin_noise_magsq))/mdl_fastDD);

  %
  % extract peaks from fastDD detection array
  %

  dfreq_fastDD = freq_fastDD(2)-freq_fastDD(1);
  nbf2 = round(f_incr_offset/dfreq_fastDD/2);
  nbf1 = 2*nbf2+1;
  nmb2 = 10;
  nmb1 = 2*nmb2+1;

  for i_chirp = 1:n_chirp

    if (df_dt_all(i_chirp)>max(df_dt_list_fastDD))
      continue;
    end

    [temp,i_freq_fastDD] = min(abs(freq_fastDD-f1_all(i_chirp)));
    [temp,i_m0] = min(abs(df_dt_list_fastDD-df_dt_all(i_chirp)));
    ii_m = max(1,i_m0-nmb2):min(size(det_fastDD,2),i_m0+nmb2);
    det_fastDD_x = det_fastDD(i_freq_fastDD+[-nbf2:nbf2],ii_m);
    freq_fastDD_x = freq_fastDD(i_freq_fastDD+[-nbf2:nbf2]);
    df_dt_fastDD_x = df_dt_list_fastDD(ii_m);
    m_list_fastDD_x = m_list_fastDD(ii_m);

    [det_max1,i1,i2] = max2d(det_fastDD_x);
    fastDD_db(i_chirp,i_snr,i_case) = 10*log10(det_max1);
    fastDD_f_est(i_chirp,i_snr,i_case) = freq_fastDD_x(i1);
    fastDD_m_est(i_chirp,i_snr,i_case) = m_list_fastDD_x(i2);
    fastDD_df_dt_est(i_chirp,i_snr,i_case) = df_dt_fastDD_x(i2);

    if isinf(snr_db)
      rel_snr_fastDD_db(i_chirp,i_case) = fastDD_db(i_chirp,i_snr,i_case) ...
                                        - bin_noise_db0 - delta_DT_fastDD_db;
    else
      SE_fastDD_db(i_chirp,i_snr,i_case) = 10*log10((det_max1 - bin_noise_magsq)/mdl_fastDD);
    end
      
    fastDD_f_error(i_chirp,i_snr,i_case) = freq_fastDD_x(i1)-f1_all(i_chirp);
    fastDD_df_dt_error(i_chirp,i_snr,i_case) = df_dt_fastDD_x(i2)-df_dt_all(i_chirp);
  end

  %
  % compute some metrics
  %

  if use_weighted_mean
    mean_taylor_db(i_snr,i_case) = GWmean(taylor_db(:,i_snr,i_case),df_dt_all);
    mean_fastDD_db(i_snr,i_case) = GWmean(fastDD_db(:,i_snr,i_case),df_dt_all);

    if isinf(snr_db)
      mean_rel_snr_taylor_db(i_N0,i_case) = GWmean(rel_snr_taylor_db(:,i_case),df_dt_all);
      mean_rel_snr_fastDD_db(i_N0,i_case) = GWmean(rel_snr_fastDD_db(:,i_case),df_dt_all);
    else
      mean_SE_taylor_db(i_snr,i_case) = GWmean(SE_taylor_db(:,i_snr,i_case),df_dt_all);
      mean_SE_fastDD_db(i_snr,i_case) = GWmean(SE_fastDD_db(:,i_snr,i_case),df_dt_all);
    end
  else
    mean_taylor_db(i_snr,i_case) = mean_nan(taylor_db(:,i_snr,i_case));
    mean_fastDD_db(i_snr,i_case) = mean_nan(fastDD_db(:,i_snr,i_case));

    if isinf(snr_db)
      mean_rel_snr_taylor_db(i_N0,i_case) = mean_nan(rel_snr_taylor_db(:,i_case));
      mean_rel_snr_fastDD_db(i_N0,i_case) = mean_nan(rel_snr_fastDD_db(:,i_case));
    else
      mean_SE_taylor_db(i_snr,i_case) = mean_nan(SE_taylor_db(:,i_snr,i_case));
      mean_SE_fastDD_db(i_snr,i_case) = mean_nan(SE_fastDD_db(:,i_snr,i_case));
    end
  end
  
  %
  % plot results as needed
  %

  if (i_case>1)
    
    plot_taylor_fastDD_output;
      
    plot_multichirp_images;
  end

end % snr
