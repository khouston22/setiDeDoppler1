%
% plot PFB bin output characteristics
%

if isinf(snr_db_list(i_snr))
  sg_limits1 = [-10 -2];
  sg_limits2 = [-8 1];
  det_limits_taylor = [-10 0];
  det_limits_fastDD = [-10 0];
  taylor_offset = 0;
  DD_offset = 0;
  z_label = 'Output dB';
  str1_fn = config_str1_fn;
  str1 = config_str1;
  str2_fn = config_str2_fn;
  str2 = config_str2;
else
  sg_limits1 = [0 10];
  sg_limits2 = [0 10];
  det_limits_taylor = [-5 5];
  det_limits_fastDD = [-5 5];
  taylor_offset = DT_taylor_db;
  DD_offset = DT_fastDD_db;
  z_label = 'Signal Excess dB';
  str1_fn = sprintf('%s-SNR-%.0f',config_str1_fn,i_snr);
  str1 = sprintf('%s, SNR=%.1f dB',config_str1,snr_db_list(i_snr));
  str2_fn = sprintf('%s-SNR-%.0f',config_str2_fn,i_snr);
  str2 = sprintf('%s, SNR=%.1f dB',config_str2,snr_db_list(i_snr));
end

if (0)
  figure(1); clf;
  %ii_freq = i_freq0 + [0:nb];
  [freq_mesh,t_mesh] = meshgrid(freq,t_out);
  plot3(freq_mesh,t_mesh,10*log10(max(1,xx)).'); % no "axis xy"
  %colormap(parula);
  %surf(freq_mesh,t_mesh,10*log10(max(1,xx)).'); % no "axis xy"
  %sg_limits1
  %colorbar; 
  zlabel('Output dB');
  xlabel('Frequency (Hz)')
  ylabel('Time (sec)')
  title(sprintf('PFB Output, %s',str1));
  %pause(2);
end

if (plot_multichirp_sg_rfi && rfi_proc_enable)
  png.file_name = sprintf('30-mc-rfi-sg-%s',str1_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  subplot(2,1,1);
  imagesc(freq,t_out,10*log10(max(1e-4,xx_rfi)).',[-4 10]); % no "axis xy"
  colorbar; colormap(parula);
  colorbar_label([],'Output dB');
  xlabel('Frequency (Hz)')
  ylabel('Time (sec)')
  title(sprintf('PFB Output, %s',str1));
  subplot(2,1,2);
  imagesc(freq,t_out,10*log10(max(1e-4,xx)).',[-4 10]); % no "axis xy"
  colorbar; colormap(parula);
  colorbar_label([],'Output dB');
  xlabel('Frequency (Hz)')
  ylabel('Time (sec)')
  title(sprintf('DD Input After RFI Processing, N-avg-t=%.0f, N-avg-f=%.0f, exp=%.0f',...
                N_rfi_avg_t,N_rfi_avg_f,rfi_weight_exp));
  png = print_fig(png);
  %pause(2);
end

if (plot_multichirp_sg)
  png.file_name = sprintf('31-mc-sg-%s',str1_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  imagesc(freq,t_out,10*log10(max(1e-4,xx)).',sg_limits1); % no "axis xy"
  colorbar; colormap(parula);
  colorbar_label([],'Output dB');
  xlabel('Frequency (Hz)')
  ylabel('Time (sec)')
  title(sprintf('PFB Output, %s',str1));
  png = print_fig(png);
  %pause(2);
end

if (plot_multichirp_cat_sg1)
  png.file_name = sprintf('32-mc-sg-%s',str1_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  imagesc(0.7:n_chirp+.3,t_out,xx_ssg_cat_db.',sg_limits2); % no "axis xy"
  %colorbar; colormap(parula);
  colorbar; colormap(jet);
  colorbar_label([],'Output dB');
  xlabel('Chirp Index')
  ylabel('Time (sec)')
  title(sprintf('DeDoppler-Shifted PFB Output, %s',str1));
  png = print_fig(png);
  %pause(2);
end

if (plot_multichirp_cat_sg2)
  png.file_name = sprintf('33-mc-sg-%s',str1_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  imagesc(0.7:n_chirp+.3,t_out,xx_ssg_cat_db.',sg_limits2); % no "axis xy"
  %colorbar; colormap(parula);
  colorbar; colormap(jet);
  colorbar_label([],'Output dB');
  if (1), v=axis; axis([0.7 min(v(2),8.5) v(3) min(v(4),100)]); end
  xlabel('Chirp Index')
  ylabel('Time (sec)')
  title(sprintf('DeDoppler-Shifted PFB Output, %s',str1));
  png = print_fig(png);
  %pause(2);
end

if plot_taylor_det_plane
  png.file_name = sprintf('36-mc-taylor-det-%s',str2_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  [temp,i_freq_t1] = min(abs(freq_taylor-(min(f1_all)-2*f_incr_offset)));
  [temp,i_freq_t2] = min(abs(freq_taylor-(max(f1_all)+2*f_incr_offset)));
  ii_f_t = i_freq_t1:i_freq_t2;
  imagesc(freq_taylor(ii_f_t),df_dt_list_taylor,10*log10(det_taylor(ii_f_t,:)).'-taylor_offset,...
               det_limits_taylor); % no "axis xy"
  colorbar; colormap(jet);
  colorbar_label([],z_label);
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Start Frequency (Hz)')
  ylabel('Drift Rate (Hz/sec)')
  title(sprintf('Taylor Multi-Chirp DeDoppler Peaks\n%s',str2));
  png = print_fig(png);
  %pause(2);
end

if plot_fastDD_det_plane
  png.file_name = sprintf('37-mc-fastDD-det-%s',str2_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  [temp,i_freq_fastDD1] = min(abs(freq_fastDD-(min(f1_all)-2*f_incr_offset)));
  [temp,i_freq_fastDD2] = min(abs(freq_fastDD-(max(f1_all)+2*f_incr_offset)));
  ii_f_fastDD = i_freq_fastDD1:i_freq_fastDD2;
  imagesc(freq_fastDD(ii_f_fastDD),df_dt_list_fastDD,10*log10(det_fastDD(ii_f_fastDD,:)).'-DD_offset,...
               det_limits_fastDD); % no "axis xy"
  colorbar; colormap(jet);
  colorbar_label([],z_label);
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Start Frequency (Hz)')
  ylabel('Drift Rate (Hz/sec)')
  title(sprintf('FastDD Multi-Chirp DeDoppler Peaks\n%s',str2));
  png = print_fig(png);
  %pause(2);
end

if 0 && plot_fastDD_det_plane
  png.file_name = sprintf('38-mc-fastDD-det-%s',str2_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  [temp,i_freq_fastDD1] = min(abs(freq_fastDD-(min(f1_all)-2*f_incr_offset)));
  [temp,i_freq_fastDD2] = min(abs(freq_fastDD-(max(f1_all)+2*f_incr_offset)));
  ii_f_fastDD = i_freq_fastDD1i_freq_fastDD2;
  [xx,yy] = meshgrid(freq_fastDD(ii_f_fastDD),df_dt_list_fastDD);
  zz = 10*log10(det_fastDD(ii_f_fastDD,:)).'-DD_offset;
  zz(find(zz<-1))=NaN;
  plot3(xx,yy,zz);
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Start Frequency (Hz)')
  ylabel('Drift Rate (Hz/sec)')
  zlabel(z_label);
  title(sprintf('FastDD Multi-Chirp DeDoppler Peaks\n%s',str2));
  png = print_fig(png);
  %pause(2);
end

if ~any(isinf(snr_db)) && plot_taylor_det_profile
  png.file_name = sprintf('34-mc-taylor-DD-det-%s',str2_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  plot(df_dt_list_fastDD,det_taylor_SE_vs_df_dt_db{i_snr,i_case},'-',...
       df_dt_list_taylor1,det_taylor_SE_vs_df_dt_db{i_snr,1},'-',...
       [df_dt_min df_dt_max],[0 0],'--k','LineWidth',0.6);
  v=axis; axis([df_dt_min df_dt_max -10 max(10*ceil(v(4)/10),10)]); v=axis;
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Drift Rate (Hz/sec)')
  ylabel('Signal Excess dB')
  legend(taylor_str0,baseline_str0,'Location','NorthEast');
  grid;
  if (EnableFinalPlots)  % plots without titles or text
    png = print_fig(png);
  end
  title(sprintf('Taylor Detection Peaks vs Baseline\n%s',str1));
  text_sc(.05,.95,sprintf('Baseline: %s: %s-%.0f',...
          run_legend_str{1},alg_baseline,N0_baseline),11);
  text_sc(.05,.05,sprintf('Mean Taylor Signal Excess: %.1f dB',...
        mean_SE_taylor_db(i_snr,i_case)),11);
  png = print_fig(png);
end

if ~any(isinf(snr_db)) && plot_fastDD_det_profile
  png.file_name = sprintf('35-mc-fastDD-DD-det-%s',str2_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  plot(df_dt_list_fastDD,det_fastDD_SE_vs_df_dt_db{i_snr,i_case},'-',...
       df_dt_list_taylor1,det_taylor_SE_vs_df_dt_db{i_snr,1},'-',...
       [df_dt_min df_dt_max],[0 0],'--k','LineWidth',2);
  v=axis; axis([df_dt_min df_dt_max -10 max(10*ceil(v(4)/10),10)]); v=axis;
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Drift Rate (Hz/sec)')
  ylabel('Signal Excess dB')
  legend(fastDD_str0,baseline_str0,'Location','NorthEast');
  grid;
  if (EnableFinalPlots)  % plots without titles or text
    png = print_fig(png);
  end
  title(sprintf('fastDD Detection Peaks vs Baseline\n%s',str1));
  text_sc(.05,.95,sprintf('Baseline: %s: %s-%.0f',...
          run_legend_str{1},alg_baseline,N0_baseline),11);
  text_sc(.05,.05,sprintf('Mean fastDD Signal Excess: %.1f dB',...
        mean_SE_fastDD_db(i_snr,i_case)),11);
  png = print_fig(png);
end

