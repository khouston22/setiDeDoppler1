%
% pfb evaluation script
%

clear;

addpath('..\misc_fns');
addpath('..\pfb_mex');

%
% parameters
%

if (1)
  %N_tap_list = [1,4,6,8,10,12];
  %N_tap_list = [1,3,4,6,8];
  N_tap_list = [1,3,4,6,8];
  %N_tap_list = [1 3];
  %N_tap_list = [3];
  window1_list = {'RectWin'};
  %window1_list = {'RectWin','Hann'};
  %window2_list = {'Hann','lpf100','lpf105','lpf110','lpf115','lpf120'};
  %window2_list = {'Hann','lpf060','lpf065','lpf070','lpf075','lpf075','lpf080','lpf085','lpf090','lpf095','lpf100'};
  window2_list = {'Hann','lpfxxx'};
  %window2_list = {'Hann'};
  M0 = 4096;
  Mx_list = [1 2];
end

use_mex = 1;
if use_mex
  fprintf('Testing PFB using mex functions\n');
  mex ../pfb_mex/pfb_Lfx_mex.c ../pfb_mex/pfb_fns.c ../pfb_mex/kiss_fft.c
end

%Lf = 2;  % enable 2:1 increase in number of bins, with overlap
fs0 = 1.0;
fs_in = M0*fs0; 
N_tap_max = max(N_tap_list);

for Mx = Mx_list
  for Lf = [1 3/2 2 4]
  %for Lf = [1 5/4 4/3 3/2 7/4 2 4]
  %for Lf = [5/4 7/4]
  %for Lf = [3/2]
  %for Lf = [4/3]
    [pLf,qLf]=rat(Lf);

    M = Mx*M0;
    if (mod(M,qLf)~=0)
      M = M*qLf/round(log2(qLf));
    end
    
    for N_tap = N_tap_list  % number of taps per subfilter

      if (N_tap==1)
        window_list = window1_list;
      else
        window_list = window2_list;
      end

      for i_window = 1:length(window_list)

        window_name = window_list{i_window};
        if strcmp(window_name(1:3),'lpf')
          window_name = choose_lpf_window(Lf,N_tap);
        end

        n_out = 1024*M0/M;
        n_in = (n_out+N_tap)*M;

        fs = fs_in/M;
        t_max = n_out/fs;

        freq = [-Lf*M/2:Lf*M/2-1]/(Lf*M)*fs_in;
        n_freq = Lf*M;
        i_freq0 = Lf*M/2+1;
        ii_freq = i_freq0 + [floor(-4*Lf):ceil(4*Lf)];

        bin_bw = fs;
        df_bin = freq(2)-freq(1);

        fprintf('M=%.0f fs=%.1f Ntap=%.0f %s %.2fx',M,fs,N_tap,window_name,Lf);

        if (1)
          f1 = -4;
          f2 =  4;
        else
          f_offset = 0;
          f1 = freq(floor(i_freq0-2*Lf))+f_offset;
          f2 = freq(ceil(i_freq0+2*Lf))+f_offset;
        end

        delta_f = f2 - f1;
        df_dt = delta_f/t_max;

        config_str = sprintf('Ntap=%.0f %s %.2fx, fs=%.1f Hz',...
          N_tap,window_name,Lf,fs);
        config_str2 = sprintf('%.2fx-Ntap-%.0f-%s-fs-%.1f',...
          Lf,N_tap,window_name,fs);

        %
        % generate chirp and run pfb
        %

        sigma = 0; A = 1;
        [x_out,freq1,t_out,coef,x_in] = gen_chirp_pfb_matrix(...
                  sigma,A,f1,df_dt,fs,n_out,M,N_tap,window_name,Lf);
        if (1)
          t_start = tic;
          [temp] = gen_chirp_pfb_matrix2(...
                   x_in,fs,n_out,M,N_tap,window_name,Lf);
          t_elapsed = toc(t_start);
          fprintf(' %.2f sec',t_elapsed*8*1024/M);
        end
        fprintf('\n');

        %
        % run PFB with noise only
        %

        sigma = sqrt(M); A = 0;
        [noise_out] = gen_chirp_pfb_matrix(...
                  sigma,A,f1,df_dt,fs,n_out,M,N_tap,window_name,Lf,use_mex);
        bin_noise_db = 20*log10(rms(noise_out(:)))';

        t_bin_ctr = (freq-f1)/df_dt;
        dt_bin_ctr = t_bin_ctr(2)-t_bin_ctr(1);

        t_bin_edge = t_bin_ctr + dt_bin_ctr/2;
        f_bin_edge = freq + df_bin/2;

        t_in = [0:n_in-1]'/fs_in;
        freq_t_out = interp1(t_in,f1+df_dt*t_in,t_out+1);

        x_out_db = 20*log10(max(1e-3,abs(x_out)));

        mean_abs_out_db = 20*log10(mean(max(abs(x_out(ii_freq,:)))));
        rms_out_db = 20*log10(rms(max(abs(x_out(ii_freq,:)))));

        snr_coh_db = mean_abs_out_db - bin_noise_db;
        snr_incoh_db = rms_out_db - bin_noise_db;

        %return;

        %
        % set up for plotting
        %

        PlotLineWidth = 1;
        AxisFontSize = 12;
        AxisFontWeight = 'bold';
        EnableFinalPlots = 1;

        output_directory=sprintf('./plots-pfb/');
        if ~exist(output_directory,'dir'), mkdir(output_directory); end

        default_figure_position = [775 43 758 591];	% set figure size for png import to pptx
        figure(1);
        set(1,'Position',default_figure_position);
        commandwindow;

        png = print_fig;
        png.output_directory = output_directory;

        plot_ID = 1;

        %
        % plot results
        %

        if (0)
          n_fft = M; min_db = -20; max_db = 10;
          [A_db,t_sg,freq_sg] = sp_gram1(x_in,0,n_fft,fs_in,'NB','NONE');
          figure(1); clf;
          imagesc(freq,t_sg,max(min_db,min(max_db,A_db))',[min_db max_db]); 
          set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
          %axis xy;
          title(sprintf('Input Spectrogram - %s',config_str));
          xlabel('Frequency (Hz)');
          ylabel('Time (sec)');
          %input('Enter any key to continue: ','s');
          pause(2);
        end

        if (0)
          png.file_name = sprintf('1-pfb-sg-%s',config_str2);
          png.file_count = plot_ID-1;
          figure(1); clf;
          imagesc(freq,t_out,x_out_db.',[-20 10]); % no "axis xy"
          set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
          colorbar; colormap(jet);
          colorbar_label([],'Output dB');
          xlabel('Frequency (Hz)')
          ylabel('Time (sec)')
          title(sprintf('PFB Output, %s',config_str));
          png = print_fig(png);
        end

        if (0)
          png.file_name = sprintf('2-pfb-ampl-%s',config_str2);
          png.file_count = plot_ID-1;
          figure(1); clf;
          plot(t_out,x_out_db(ii_freq,:),'-'); hold on;
          v=axis; axis([0 t_max -30 10]); v=axis; 
          %plot([1 1]'*t_bin_ctr(ii_freq),v(3:4)'*ones(size(ii_freq)),'--'); hold on;
          plot([1 1]'*t_bin_edge(ii_freq),v(3:4)'*ones(size(ii_freq)),'-.','LineWidth',PlotLineWidth); 
          set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
          xlabel('Time (sec)')
          ylabel('Amplitude (dB)')
          title(sprintf('PFB Time Series, %s',config_str));
          text_sc(.05,.95,sprintf('Envelope Mean=%.2f dB, RMS=%.2f dB',mean_abs_out_db,rms_out_db));
          text_sc(.05,.90,sprintf('Bin Noise=%.2f dB',bin_noise_db));
          %text_sc(.05,.85,sprintf('Relative Coherent SNR=%.2f dB, Incoherent SNR=%.2f dB',snr_coh_db,snr_incoh_db));
          text_sc(.05,.85,sprintf('Relative SNR=%.2f dB',snr_incoh_db));
          grid;
          png = print_fig(png);
        end

        if (1)
          png.file_name = sprintf('3-pfb-ampl-%s',config_str2);
          png.file_count = plot_ID-1;
          figure(1); clf;
          plot(freq_t_out,x_out_db(ii_freq,:),'-','LineWidth',PlotLineWidth); hold on;
          v=axis; axis([-1.5 1.5 -30 10]); v=axis; 
          plot([1 1]'*f_bin_edge(ii_freq),v(3:4)'*ones(size(ii_freq)),'-.','LineWidth',PlotLineWidth); 
          set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
          xlabel('Frequency (Hz)')
          ylabel('Amplitude (dB)')
          title(sprintf('PFB Bin Freq Response, %s',config_str));
          grid;
          %png = print_fig(png);
          text_sc(.05,.95,sprintf('Envelope Mean=%.2f dB, RMS=%.2f dB',mean_abs_out_db,rms_out_db));
          text_sc(.05,.90,sprintf('Bin Noise=%.2f dB',bin_noise_db));
          %text_sc(.05,.85,sprintf('Relative Coherent SNR=%.2f dB, Incoherent SNR=%.2f dB',snr_coh_db,snr_incoh_db));
          text_sc(.05,.85,sprintf('Relative SNR=%.2f dB',snr_incoh_db));
          png = print_fig(png);
        end

        if (0)
          png.file_name = sprintf('4-pfb-phase-%s',config_str2);
          png.file_count = plot_ID-1;
          figure(1); clf;
          plot(freq_t_out,angle(x_out(ii_freq,:))*180/pi + 400*(ii_freq-i_freq0)','-'); hold on;
          v=axis; axis([f1 f2 v(3:4)]); v=axis; 
          set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
          plot([1 1]'*f_bin_edge(ii_freq),v(3:4)'*ones(size(ii_freq)),'-.','LineWidth',PlotLineWidth); 
          xlabel('Frequency (Hz)')
          ylabel('Phase (deg)')
          title(sprintf('PFB Phase Time Series, %s',config_str));
          grid;
          png = print_fig(png);
        end

        if (0)
          % examine phase compensation for SCP detection
          png.file_name = sprintf('5-pfb-conj-product-%s',config_str2);
          png.file_count = plot_ID-1;
          figure(1); clf;
          x1 = x_out(i_freq0+[-1:1],:);
          w1 = exp(-1j*4*pi*([0:n_out-1]'.^2)/n_out).';
          %w1 = exp(-1j*4*pi*(([0:n_out-1]'-n_out/2).^2)/n_out).';
          x1 = x1.*(ones(3,1)*w1);
          ii1d = [1:n_out];
          ii1a = [2:n_out n_out];
          %ii1d = [3:n_out n_out n_out];
          cp1 = x1(:,ii1a).*conj(x1(:,ii1d));
          subplot(2,1,1); 
          plot(freq_t_out,abs(cp1),'LineWidth',PlotLineWidth); 
          set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
          xlabel('Frequency (Hz)')
          ylabel('Amplitude')
          title(sprintf('PFB Conj Product Amplitude, %s',config_str));
          grid; 
          subplot(2,1,2); 
          plot(freq_t_out,angle(cp1)*180/pi,'LineWidth',PlotLineWidth); 
          set(gca,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
          grid; 
          xlabel('Frequency (Hz)')
          ylabel('Phase (deg)')
          title(sprintf('PFB Conj Product Phase'));
          png = print_fig(png);
        end

      end % window
    end % N_tap
  end % Lf

end % MX
