%
% initialize arrays for a run
%

taylor_best_db = NaN(n_chirp,n_snr,n_case);
bin_noise_db = NaN(n_chirp,n_snr,n_case);
rel_snr_best_db = NaN(n_chirp,n_case);

% scp_db = NaN(n_chirp,max_n_diff,n_case);
% mean_scp_db = NaN(max_n_diff,n_case);

taylor_db = NaN(n_chirp,n_snr,n_case);
taylor_SE_db = NaN(n_chirp,n_snr,n_case);
taylor_m_est = NaN(n_chirp,n_snr,n_case);
taylor_df_dt_est = NaN(n_chirp,n_snr,n_case);
taylor_f_est = NaN(n_chirp,n_snr,n_case);
rel_snr_taylor_db = NaN(n_chirp,n_case);
SE_taylor_db = NaN(n_chirp,n_snr,n_case);

fastDD_db = NaN(n_chirp,n_snr,n_case);
fastDD_SE_db = NaN(n_chirp,n_snr,n_case);
fastDD_m_est = NaN(n_chirp,n_snr,n_case);
fastDD_df_dt_est = NaN(n_chirp,n_snr,n_case);
fastDD_f_est = NaN(n_chirp,n_snr,n_case);
rel_snr_fastDD_db = NaN(n_chirp,n_case);
SE_fastDD_db = NaN(n_chirp,n_snr,n_case);

baseline_db = NaN(n_chirp,n_snr);
baseline_SE_db = NaN(n_chirp,n_snr,n_case);
rel_snr_baseline_db = NaN(n_chirp,1);
SE_baseline_db = NaN(n_chirp,n_snr);

mean_best_db = NaN(1,n_case);
mean_rel_snr_best_db = NaN(1,n_case);

mean_taylor_db = NaN(n_snr,n_case);
mean_SE_taylor_db = NaN(n_snr,n_case);

mean_fastDD_db = NaN(n_snr,n_case);
mean_SE_fastDD_db = NaN(n_snr,n_case);

mean_baseline_db = NaN(1,1);
mean_baseline_SE_db = NaN(1,n_snr);
mean_rel_snr_baseline_db = NaN(1,1);
mean_SE_baseline_db = NaN(1,n_snr);

t_cpu_taylor = NaN(n_chirp,n_snr,n_case);
t_cpu_fastDD = NaN(n_chirp,n_snr,n_case);

if (i_N0==1)
  delta_snr_taylor_db = NaN(n_N0,n_case);
  delta_snr_fastDD_db = NaN(n_N0,n_case);
  
  mean_rel_snr_taylor_db = NaN(n_N0,n_case);
  mean_rel_snr_fastDD_db = NaN(n_N0,n_case);

  mean_t_cpu_pfb = NaN(n_N0,n_case);
  mean_t_cpu_taylor = NaN(n_N0,n_case);
  mean_t_cpu_fastDD = NaN(n_N0,n_case);
  mean_t_cpu_taylor_total = NaN(n_N0,n_case);
  mean_t_cpu_fastDD_total = NaN(n_N0,n_case);

  det_snr_taylor_db = NaN(n_N0,n_case);
  det_snr_fastDD_db = NaN(n_N0,n_case);
  det_snr_baseline_db = NaN(n_N0,n_case);
  delta_det_snr_taylor_db = NaN(n_N0,n_case);
  delta_det_snr_fastDD_db = NaN(n_N0,n_case);
  
  value_metric_taylor = NaN(n_N0,n_case);
  value_metric_fastDD = NaN(n_N0,n_case);
end




