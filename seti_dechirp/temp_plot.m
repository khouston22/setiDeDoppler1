  figure(1); clf;
  subplot(2,1,1);
  imagesc(freq,t_out,10*log10(max(1e-4,xx)).',[-4 10]); % no "axis xy"
  colorbar; colormap(parula);
  colorbar_label([],'Output dB');
  xlabel('Frequency (Hz)')
  ylabel('Time (sec)')
  title(sprintf('PFB Output, %s',str1));
  subplot(2,1,2);
  imagesc(freq,t_out,10*log10(max(1e-4,xx_agc)).',[-4 10]); % no "axis xy"
  colorbar; colormap(parula);
  colorbar_label([],'Output dB');
  xlabel('Frequency (Hz)')
  ylabel('Time (sec)')
  title(sprintf('DD Input, %s',str1));
