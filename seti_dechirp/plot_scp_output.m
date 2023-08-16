%
% plot scp outputs
%

if (0)
  png.file_name = sprintf('21-wdx-%s',sig_config_str2_fn);
  png.file_count = plot_ID-1;
  figure(1); clf;
  for n_diff = [1:max_n_diff]
    subplot(max_n_diff,2,2*n_diff-1)
    plot(t_out,10*log10(abs(wdx_maxi(:,n_diff))));
    v=axis; axis([0 t_max -4 2]); v=axis;
    xlabel('Time (sec)')
    ylabel('Amplitude dB')
    if n_diff==1
      title(sprintf('%s',config_str2));
    end
    text_sc(.05,.9,sprintf('n-diff=%.0f scp=%.2f dB',n_diff,scp_db(n_diff)));
    grid;
    subplot(max_n_diff,2,2*n_diff)
    plot(t_out,angle(wdx_maxi(:,n_diff))*180/pi);
    v=axis; axis([0 t_max -200 200]); v=axis;
    xlabel('Time (sec)')
    ylabel('Degrees')
    text_sc(.05,.9,sprintf('n-diff=%.0f',n_diff));
    grid;
  end
  png = print_fig(png);
end
      
