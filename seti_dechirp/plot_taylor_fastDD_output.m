%
% plot taylor outputs
%

if ~(plot_taylor_det_peaks||plot_taylor_shifted_spgram||...
     plot_fastDD_det_peaks||plot_fastDD_shifted_spgram)
  return;
end

n_plot = min(3,n_chirp);
ii_chirp = 1:round(n_chirp/n_plot):n_chirp;
for i_chirp = ii_chirp
  
  f1 = f1_all(i_chirp);
  f2 = f2_all(i_chirp);
  m0_chirp = m0_chirp_all(i_chirp);
  df_dt = df_dt_all(i_chirp);

  sig_str = sprintf('Chirp %.2f..%.2f Hz (%.1f bins in %.0f sec, %.3f Hz/sec)',...
    f1,f2,m0_chirp,t_max,df_dt);
  sig_str_fn = sprintf('f_%.2f-%.1fBins-%.2fHzsec-%.1fsec',...
    f1,m0_chirp,t_max,df_dt);
  if ~isempty(mod_str)
    sig_str = sprintf('%s\n%s',sig_str,mod_str);
    sig_str_fn = [sig_str_fn '-' mod_str_fn];
  end

  sig_config_str0 = sprintf('%s\n%s',sig_str,config_str0);
  sig_config_str1 = sprintf('%s\n%s',sig_str,config_str1);
  sig_config_str2 = sprintf('%s\n%s',sig_str,config_str2);

  sig_config_str0_fn = sprintf('%s-%s',sig_str_fn,config_str0_fn);
  sig_config_str1_fn = sprintf('%s-%s',sig_str_fn,config_str1_fn);
  sig_config_str2_fn = sprintf('%s-%s',sig_str_fn,config_str2_fn);

  taylor_string1 = sprintf('Taylor1 N0=%.0f freq=%.2f m=%4.0f dfdt=%.3f A=%6.2f dB relSNR=%6.2f dB',...
    N0,taylor_f_est(i_chirp,i_snr,i_case),taylor_m_est(i_chirp,i_snr,i_case),taylor_df_dt_est(i_chirp,i_snr,i_case),...
    taylor_db(i_chirp,i_snr,i_case),rel_snr_taylor_db(i_chirp,i_case));
  fprintf('%s\n',taylor_string1);
  
  fastDD_string1 = sprintf('fastDD N0=%.0f freq=%.2f m=%4.0f dfdt=%.3f A=%6.2f dB relSNR=%6.2f dB',...
    N0,fastDD_f_est(i_chirp,i_snr,i_case),fastDD_m_est(i_chirp,i_snr,i_case),fastDD_df_dt_est(i_chirp,i_snr,i_case),...
    fastDD_db(i_chirp,i_snr,i_case),rel_snr_fastDD_db(i_chirp,i_case));
  fprintf('%s\n',fastDD_string1);

  if (plot_fastDD_det_peaks)
    png.file_name = sprintf('21-det-output-%s-fastDD',sig_config_str2_fn);
    png.file_count = plot_ID-1;
    for i_plot = 1:2
      if (i_plot==2 && ~EnableFinalPlots)
        break;
      end
      figure(1); clf;
      subplot(2,1,1);
      [temp,i_freq] = min(abs(f1-freq_fastDD));
      nf = length(freq_fastDD);
      ii_freq1 = max(1,i_freq-Lr*4):min(nf,i_freq+Lr*4);
      ii_freq2 = max(1,i_freq-Lr*4):min(nf,i_freq+Lr*4);
      if (pos_neg==1)
        ii_m = max(1,1 + round(Lr*m0_chirp) -Lr*10):min(Lr*Nt_pfb,1 + round(Lr*m0_chirp)+Lr*10);
      else
        ii_m = max(1,Lr*Nt_pfb + 1 + round(Lr*m0_chirp) -Lr*10):min(Lr*Nt_pfb,Lr*Nt_pfb + 1 + round(Lr*m0_chirp)+Lr*10);
      end
      imagesc(freq_fastDD(ii_freq1),df_dt_list_fastDD(ii_m),10*log10(det_fastDD(ii_freq1,ii_m)).',...
                   [-20 0]); % no "axis xy"
      colorbar; colormap(jet);
      colorbar_label([],'Output dB');
      set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
      xlabel('Start Frequency (Hz)')
      ylabel('Drift Rate (Hz/sec)')
      if (i_plot==1)
        title(sprintf('FastDD DeDoppler Peak, %s',sig_config_str2));
      else
        title(sprintf('FastDD DeDoppler Peak, %.2fx Lr=%.2f',Lf,Lr));
      end
      subplot(2,2,3);
      plot(freq_fastDD(ii_freq1),10*log10(det_fastDD(ii_freq1,ii_m)),'-*'); 
      v=axis; axis([v(1:2) -10 2]);
      set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
      xlabel('Start Frequency (Hz)')
      ylabel('Amplitude dB');
      title('Frequency Profile');
      if (i_plot==1)
        text_sc(.05,.95,sprintf('Peak=%.2f dB vs %.2f dB Ideal',...
                fastDD_db(i_chirp,i_snr,i_case),taylor_best_db(i_chirp,i_snr,i_case)));
        text_sc(.05,.88,sprintf('f1=%.2f Est vs %.2f Hz',...
                fastDD_f_est(i_chirp,i_snr,i_case),f1));
      end
      grid;
      subplot(2,2,4);
      plot(df_dt_list_fastDD(ii_m),10*log10(det_fastDD(ii_freq2,ii_m)),'-*'); 
      v=axis; axis([v(1:2) -10 2]);
      set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
      xlabel('Drift Rate (Hz/sec)')
      ylabel('Amplitude dB');
      title('Drift Rate Profile');
      if (i_plot==1)
        text_sc(.05,.95,sprintf('Execution %.1f sec',t_cpu_fastDD(i_chirp,i_snr,i_case)));
        text_sc(.05,.88,sprintf('dF/dt=%.3f Est vs %.3f Hz/sec',...
                fastDD_df_dt_est(i_chirp,i_snr,i_case),df_dt));
      end
      grid;
      png = print_fig(png);
    end
  end      

  if (plot_taylor_det_peaks)
    png.file_name = sprintf('22-det-output-%s-taylor',sig_config_str2_fn);
    png.file_count = plot_ID-1;
    for i_plot = 1:2
      if (i_plot==2 && ~EnableFinalPlots)
        break;
      end
      figure(1); clf;
      subplot(2,1,1);
      [temp,i_freq] = min(abs(f1-freq_taylor));
      nf = length(freq_taylor);
      ii_freq1 = max(1,i_freq-6):min(nf,i_freq+6);
      ii_freq2 = max(1,i_freq-4):min(nf,i_freq+4);
      if (pos_neg==1)
        ii_m = 1 + round(m0_chirp) + [-10:10];
      else
        ii_m = Nt_pfb + 1 + round(m0_chirp) + [-10:10];
      end
      imagesc(freq_taylor(ii_freq1),df_dt_list_taylor(ii_m),...
                   10*log10(det_taylor(ii_freq1,ii_m)).',[-20 0]); % no "axis xy"
      colorbar; colormap(jet);
      colorbar_label([],'Output dB');
      set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
      xlabel('Start Frequency (Hz)')
      ylabel('Drift Rate (Hz/sec)')
      if (i_plot==1)
        title(sprintf('Taylor DeDoppler Peak, %s',sig_config_str2));
      else
        title(sprintf('Taylor DeDoppler Peak, %.2fx Lr=%.2f',Lf,Lr));
      end
      subplot(2,2,3);
      plot(freq_taylor(ii_freq1),10*log10(det_taylor(ii_freq1,ii_m)),'-*'); 
      v=axis; axis([v(1:2) -10 2]);
      set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
      xlabel('Start Frequency (Hz)')
      ylabel('Amplitude dB');
      title('Frequency Profile');
      if (i_plot==1)
        text_sc(.05,.95,sprintf('Peak=%.2f dB vs %.2f dB Ideal',...
                taylor_db(i_chirp,i_snr,i_case),taylor_best_db(i_chirp,i_snr,i_case)));
        text_sc(.05,.88,sprintf('f1=%.2f Est vs %.2f Hz',...
                taylor_f_est(i_chirp,i_snr,i_case),f1));
      end
      grid;
      subplot(2,2,4);
      plot(df_dt_list_taylor(ii_m),10*log10(det_taylor(ii_freq2,ii_m)),'-*'); 
      v=axis; axis([v(1:2) -10 2]);
      set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
      xlabel('Drift Rate (Hz/sec)')
      ylabel('Amplitude dB');
      title('Drift Rate Profile');
      if (i_plot==1)
        text_sc(.05,.95,sprintf('Execution %.1f sec',t_cpu_taylor(i_chirp,i_snr,i_case)));
        text_sc(.05,.88,sprintf('dF/dt=%.3f Est vs %.3f Hz/sec',...
                taylor_df_dt_est(i_chirp,i_snr,i_case),df_dt));
      end
      grid;
      png = print_fig(png);
    end      
  end
  
  if (0)
    png.file_name = sprintf('26-taylorDD-nstage-%.0f-%s',n_stage,sig_config_str2_fn);
    png.file_count = plot_ID-1;
    figure(1); clf;
    subplot(2,1,1);
    [temp1,i_mc] = min(abs(m_list-mc));
    ii_mc = max(1,i_mc-15):min(i_mc+15,length(m_list));
    if (0)
      imagesc(freq(ii_freq2),t_out,x_out_db(ii_freq2,:).',[-20 10]); % no "axis xy"
      colorbar; colormap(jet);
      colorbar_label([],'Output dB');
      xlabel('Frequency (Hz)')
      ylabel('Time (sec)')
    else
      ii_mc2 = max(1,i_mc-5):min(i_mc+5,length(m_list));
      ii_freq4 = i_freq0 + round((f1+f2)/2/df_bin) +[-4:4];
      plot(m_list(ii_mc2),10*log10(det_taylor(ii_freq4,ii_mc2)),'-*');
      v=axis; axis([v(1:2) -10 2]);
      grid;
      xlabel('Drift Index')
      ylabel('dB')
      title(sprintf('Energy Peak, %s',sig_config_str2));
    end
    subplot(2,1,2);
    ii_freq3 = i_freq0 + round((f1+f2)/2/df_bin) +[-25:25];
    imagesc(freq(ii_freq3),m_list(ii_mc),10*log10(det_taylor(ii_freq3,ii_mc)).',[-20 10]); % no "axis xy"
    colorbar; colormap(jet);
    colorbar_label([],'TaylorDD dB');
    xlabel('Frequency (Hz)')
    ylabel('Drift index')
    title(sprintf('TaylorDD Output: %s',taylor_string1));
    png = print_fig(png);
    pause(2);
  end
  
end

