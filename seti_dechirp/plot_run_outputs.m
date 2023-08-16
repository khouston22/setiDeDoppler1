%
% plot final taylor and fastDD outputs
%

snr_db = snr_db_list(i_snr);

if (0 & isinf(snr_db))
  if isempty(mod_str)
    png.file_name = sprintf('06-relSNR-db-%s-taylor',config_str2_fn);
  else
    png.file_name = sprintf('06-relSNR-db-%s-%s-taylor',config_str2_fn,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
%   plot(df_dt_all*Ts0.^2,rel_snr_best_db(:,i_case),'-*',...
%        df_dt_all*Ts0.^2,rel_snr_taylor_db(:,i_case),'-^',...
%        df_dt_all*Ts0.^2,rel_snr_baseline_db(:),'-o','LineWidth',PlotLineWidth);
  plot(df_dt_all,rel_snr_best_db(:,i_case),'-*',...
       df_dt_all,rel_snr_taylor_db(:,i_case),'-^',...
       df_dt_all,rel_snr_baseline_db(:),'-o','LineWidth',PlotLineWidth);
  v=axis; axis([v(1:2) -3 4]); v=axis;
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Drift Rate (Hz/sec)');
  ylabel('Relative SNR dB');
  legend('Best Possible',taylor_str0,'Baseline');
  grid;
  if (EnableFinalPlots)  % plots without titles or text
    png = print_fig(png);
  end
  title(sprintf('Taylor DeDoppler Relative SNR, %s',config_str1));
  text_sc(.05,.95,best_str);
  text_sc(.05,.92,taylor_str1);
  text_sc(.05,.89,baseline_str1);
  text_sc(.05,.86,t_cpu_str);
  text_sc(.05,.83,noise_str);
  text_sc(.05,.80,DT_str);
  if ~isempty(mod_str)
    text_sc(.05,.77,mod_str);
  end
  png = print_fig(png);
end

if (1 & isinf(snr_db))
  if isempty(mod_str)
    png.file_name = sprintf('06-relSNR-db-%s-fastDD',config_str2_fn);
  else
    png.file_name = sprintf('06-relSNR-db-%s-%s-fastDD',config_str2_fn,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  plot(df_dt_all,rel_snr_best_db(:,i_case),'-*',...
       df_dt_all,rel_snr_fastDD_db(:,i_case),'-^',...
       df_dt_all,rel_snr_baseline_db(:),'-o','LineWidth',PlotLineWidth);
  v=axis; axis([v(1:2) -3 4]); v=axis;
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Drift Rate (Hz/sec)');
  ylabel('Relative SNR dB');
  legend('Best Possible',fastDD_str0,'Baseline');
  grid;
  if (EnableFinalPlots)  % plots without titles or text
    png = print_fig(png);
  end
  title(sprintf('FastDD DeDoppler Relative SNR, %s',config_str1));
  text_sc(.05,.95,best_str);
  text_sc(.05,.92,fastDD_str1);
  text_sc(.05,.89,baseline_str1);
  text_sc(.05,.86,t_cpu_str);
  text_sc(.05,.83,noise_str);
  text_sc(.05,.80,DT_str);
  if ~isempty(mod_str)
    text_sc(.05,.77,mod_str);
  end
  png = print_fig(png);
end

if (1 & isinf(snr_db))
  if isempty(mod_str)
    png.file_name = sprintf('07-relSNR-db-%s',config_str2_fn);
  else
    png.file_name = sprintf('07-relSNR-db-%s-%s',config_str2_fn,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  plot(df_dt_all,rel_snr_best_db(:,i_case),'-*',...
       df_dt_all,rel_snr_fastDD_db(:,i_case),'-^',...
       df_dt_all,rel_snr_taylor_db(:,i_case),'-^',...
       df_dt_all,rel_snr_baseline_db(:),'-o','LineWidth',PlotLineWidth);
  v=axis; axis([v(1:2) -3 4]); v=axis;
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Drift Rate (Hz/sec)');
  ylabel('Relative SNR dB');
  legend('Best Possible',fastDD_str0,taylor_str0,'Baseline')
  grid;
  if (EnableFinalPlots)  % plots without titles or text
    png = print_fig(png);
  end
  title(sprintf('FastDD & Taylor DeDoppler Relative SNR, %s',config_str1));
  text_sc(.05,.95,best_str);
  text_sc(.05,.92,fastDD_str1);
  text_sc(.05,.89,taylor_str1);
  text_sc(.05,.86,baseline_str1);
  text_sc(.05,.83,t_cpu_str);
  text_sc(.05,.80,noise_str);
  text_sc(.05,.77,DT_str);
  if ~isempty(mod_str)
    text_sc(.05,.74,mod_str);
  end
  png = print_fig(png);
end

if (0 & isinf(snr_db))
  if isempty(mod_str)
    png.file_name = sprintf('08-relSNR-db-%s',config_str2_fn);
  else
    png.file_name = sprintf('08-relSNR-db-%s-%s',config_str2_fn,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  plot(mod_ctr(m0_chirp_all,1),rel_snr_best_db(:,i_case),'*',...
       mod_ctr(m0_chirp_all,1),rel_snr_fastDD_db(:,i_case),'^','LineWidth',PlotLineWidth);
  v=axis; axis([v(1:2) -3 4]); v=axis;
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Drift Rate Fractional Index');
  ylabel('Relative SNR dB');
  legend('Best Possible','fastDD')
  grid;
  if (EnableFinalPlots)  % plots without titles or text
    png = print_fig(png);
  end
  title(sprintf('FastDD DeDoppler Relative SNR, %s',config_str1));
  text_sc(.05,.95,best_str);
  text_sc(.05,.92,fastDD_str1);
  text_sc(.05,.89,baseline_str1);
  text_sc(.05,.86,noise_str);
  text_sc(.05,.83,DT_str);
  if ~isempty(mod_str)
    text_sc(.05,.80,mod_str);
  end
  png = print_fig(png);
end


if (0 & ~isinf(snr_db))
  if isempty(mod_str)
    png.file_name = sprintf('13-SE-db-%s-SNR-%.0f',config_str1_fn,i_snr);
  else
    png.file_name = sprintf('13-SE-db-%s-%s-SNR-%.0f',config_str1_fn,mod_str_fn,i_snr);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  plot(df_dt_all,SE_fastDD_db(:,i_snr,i_case),'-^',...
       df_dt_all,SE_taylor_db(:,i_snr,i_case),'-o',...
       df_dt_all,SE_baseline_db(:,i_snr),'-v','LineWidth',PlotLineWidth);
  v=axis; axis([v(1:2) -10 10]); v=axis;
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Drift Rate (Hz/sec)');
  ylabel('Signal Excess dB');
  legend(fastDD_str0,taylor_str0,baseline_str0,'Location','NorthEast');
  grid;
  if (EnableFinalPlots)  % plots without titles or text
    png = print_fig(png);
  end
  title(sprintf('Taylor/FastDD DeDoppler Signal Excess, %s, SNR=%.1f dB',config_str1,snr_db));
  text_sc(.05,.95,taylor_str2);
  text_sc(.05,.92,fastDD_str2);
  text_sc(.05,.89,baseline_str2);
  if ~isempty(mod_str)
    text_sc(.05,.86,mod_str);
  end
  png = print_fig(png);
end

if (1 & ~any(isinf(snr_db)) & (i_snr==n_snr))
  if isempty(mod_str)
    png.file_name = sprintf('14-taylor-fastDD-det-db-%s',config_str1_fn);
  else
    png.file_name = sprintf('14-taylor-fastDD-det-db-%s-%s',config_str1_fn,mod_str_fn);
  end
  png.file_count = plot_ID-1;
  figure(1); clf;
  plot(snr_db_list,mean_SE_fastDD_db(:,i_case),'-^',...
       snr_db_list,mean_SE_taylor_db(:,i_case),'-v',...
       snr_db_list,mean_SE_baseline_db,'-o',...
       [min(snr_db_list) max(snr_db_list)],[0 0],'--k','LineWidth',PlotLineWidth);
  v=axis; axis([v(1:2) -4 max(snr_db_list(end)+2,v(4))]); v=axis;
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Baseline Input SNR dB');
  ylabel('Signal Excess dB');
  legend(fastDD_str0,taylor_str0,baseline_str0,'Location','SouthEast');
  grid;
  if (EnableFinalPlots)  % plots without titles or text
    png = print_fig(png);
  end
  title(sprintf('Taylor/fastDD Signal Excess vs SNR, %s',config_str1));
  text_sc(.05,.95,sprintf('fastDD-%.0f %s: Det SNR = %.2f dB, %6.2f dB Over Baseline',...
    N0,run_legend_str{i_case},det_snr_fastDD_db(i_N0,i_case),...
    delta_det_snr_fastDD_db(i_N0,i_case)));
  text_sc(.05,.92,sprintf('Taylor-%.0f   %s: Det SNR = %.2f dB, %6.2f dB Over Baseline',...
    N0_taylor,run_legend_str{i_case},det_snr_taylor_db(i_N0,i_case),...
    delta_det_snr_taylor_db(i_N0,i_case)));
  text_sc(.05,.89,sprintf('Baseline    %s: Det SNR = %.2f dB (%s-%.0f)',...
    run_legend_str{1},det_snr_baseline_db,alg_baseline,N0_baseline));
  print_fig(png);
end

