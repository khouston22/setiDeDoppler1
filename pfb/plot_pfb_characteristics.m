%
% plot PFB bin output characteristics
%

PlotLineWidth = .6;
AxisFontSize = 12;
AxisFontWeight = 'bold';
EnableFinalPlots = 1;

if (plot_pfb_sg)
  png.file_name = sprintf('11-pfb-sg-%s',sig_config_str1_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  ii_freq = i_freq0 + [0:nb];
  if isinf(snr_db)
    imagesc(freq(ii_freq),t_out,x_out_db(ii_freq,:).',[-10 5]); % no "axis xy"
    colorbar; colormap(jet);
  else
    imagesc(freq(ii_freq),t_out,x_out_db(ii_freq,:).',[0 10]); % no "axis xy"
    colorbar; colormap(parula);
  end
  v=axis; axis([0 df_dt*Ts*50*1.1 t_out([1 50])])
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  colorbar_label([],'Output dB');
  xlabel('Frequency (Hz)')
  ylabel('Time (sec)')
  title(sprintf('PFB Output, %s',sig_config_str1));
  png = print_fig(png);
  %pause(2);
end

if (0)
  png.file_name = sprintf('14-pfb-3D-%s',sig_config_str1_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  ii_freq = i_freq0 + [0:nb];
  if isinf(snr_db)
    plot3(freq(ii_freq)'*ones(1,n_out),ones(nb+1,1)*t_out,x_out_db(ii_freq,:)); % no "axis xy"
    v=axis; axis([v(1:4) -5 5]);
  else
    plot3(freq(ii_freq)'*ones(1,n_out),ones(nb+1,1)*t_out,x_out_db(ii_freq,:)); % no "axis xy"
    v=axis; axis([v(1:4) 0 10]);
  end
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Frequency (Hz)')
  ylabel('Time (sec)')
  zlabel('Output dB');
  title(sprintf('PFB Output, %s',sig_config_str1));
  png = print_fig(png);
  %pause(2);
end

if (0)
  png.file_name = sprintf('12-adjacent-pfb-ampl-%s',sig_config_str1_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  ii_freq = i_freq0 + [0:nb2];
  plot(freq(ii_freq),xx(ii_freq,1:20:n_out/2),'-*'); hold on;
  plot([min(freq(ii_freq)) max(freq(ii_freq))],[1 1],'--k','LineWidth',1);
  v=axis; axis([min(freq(ii_freq)) max(freq(ii_freq)) 0 1.8]); 
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Frequency (Hz)')
  ylabel('Mag Square Amplitude')
  title(sprintf('PFB Adjacent Bin Amplitudes, %s',sig_config_str1));
  grid;
  text_sc(.02,.95,'Bin Amplitudes Shown for Every 20 Time Points');
  text_sc(.02,.92,sprintf('Best mean = %.2f dB',turbo_best_db));
  text_sc(.02,.89,['Bin Noise = ' sprintf('%.2f ',bin_noise_db) 'dB']);
  text_sc(.02,.86,['Relative SNR = ' sprintf('%.2f ',rel_snr_best_db) 'dB']);
  png = print_fig(png);
end

if (plot_pfb_peaks)
  png.file_name = sprintf('13-peak-pfb-ampl-%s',sig_config_str1_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  n_max = floor(n_out/4);
  plot(t_out(1:n_max),maxk(xx(:,1:n_max),3),'-*','LineWidth',PlotLineWidth); hold on;
  plot(t_out([1 n_max]),[1 1],'--k','LineWidth',1);
  v=axis; axis([t_out([1 n_max]) 0 1.8]); 
  set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
  xlabel('Time (sec)')
  ylabel('Max Bin Amplitude (Mag. Squared)')
  title(sprintf('PFB Peak Bin Amplitudes, %s',sig_config_str1));
  text_sc(.02,.95,sprintf('Best mean = %.2f dB',turbo_best_db));
  text_sc(.02,.92,['Bin Noise = ' sprintf('%.2f ',bin_noise_db) 'dB']);
  text_sc(.02,.89,['Relative SNR = ' sprintf('%.2f ',rel_snr_best_db) 'dB']);
  grid;
  legend('Maximum','2nd Highest','3rd Highest','Nominal','Location','NorthEast');
  png = print_fig(png);
end

