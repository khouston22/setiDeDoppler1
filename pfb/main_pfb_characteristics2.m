%
% pfb evaluation script
%
cd
clear;

addpath('..\misc_fns');
addpath('..\pfb_mex');

%
% parameters
%

use_mex = 1;
if use_mex
  fprintf('Testing PFB using mex functions\n');
  mex ../pfb_mex/pfb_Lfx_mex.c ../pfb_mex/pfb_fns.c ../pfb_mex/kiss_fft.c
end

%N_tap_list = [1,4,6,8];
N_tap_list = [1,3];

Lf_list = [1 3/2 2 4];

window1_list = {'RectWin'};
%window1_list = {'RectWin','Hann'};
window2_list = {'Hann','lpfxxx'};

%m_chirp_list = reshape([10 100 250 500 750 1000]+[0 .5]',1,[]);
%m_chirp_list = [100 250 500 750]+.25;
m_chirp_list = [100 400]+.25;
f_offset_list = 0;

t_max0 = 1024;
fs0 = 1;
Mx_list = [1 2];

i_chirp = 0;
for m_chirp = m_chirp_list
  for f_offset = f_offset_list
    i_chirp = i_chirp + 1;
    m_chirp_all(i_chirp) = m_chirp;
    f_offset_all(i_chirp) = f_offset;
    df_dt_all(i_chirp) = m_chirp*fs0/t_max0;
  end
end
n_chirp = length(m_chirp_all);

n_case = length(Lf_list)*length(find(N_tap_list==1))*length(window1_list);
n_case = n_case + length(Lf_list)*length(find(N_tap_list>1))*length(window2_list);
n_case = n_case*length(Mx_list);
i_case = 0;

mod_type = []; n_mod_ID = 1;
%mod_type='BPSK'; n_mod_ID = 10;
%mod_type='AM';   n_mod_ID = 10;

%
% set up for plotting
%

output_directory=sprintf('./plots-pfb-env2/');
if ~exist(output_directory,'dir'), mkdir(output_directory); end

default_figure_position = [775 43 758 591];	% set figure size for png import to pptx
figure(1);
set(1,'Position',default_figure_position);
commandwindow;

png = print_fig;
png.output_directory = output_directory;

%
% nested loops: Modulation, Lf, N_tap
%

for mod_ID=1:n_mod_ID
  mod_params = [];
  if isempty(mod_type)
    mod_str = '';
  else
    mod_params.mod_type = mod_type;
    mod_params.f_sym = (.25*(mod_ID-1)+.01)*fs;
    mod_params.T_sym=1/mod_params.f_sym;
    if strcmp(mod_params.mod_type,'BPSK')
      mod_str = sprintf('BPSK, f-sym=%.2f, T-sym=%.2f',mod_params.f_sym,...
                         mod_params.T_sym);
      mod_str_fn = sprintf('BPSK-%.0f',mod_ID);
    elseif strcmp(mod_params.mod_type,'AM')
      mod_params.Am=.25;
      mod_str = sprintf('AM %.2f, f-sym=%.2f, T-sym=%.2f',mod_params.Am,...
                         mod_params.f_sym,mod_params.T_sym);
      mod_str_fn = sprintf('AM-%.0f',mod_ID);
    end
    if (mod_ID==1)
      mod_str = 'No Modulation';
      mod_params.mod_type = [];
    end
  end
  
  plot_ID = mod_ID;
  
  for Lf = Lf_list

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

        [pLf,qLf]=rat(Lf);
        if mod(qLf,2)==0
          M0 = 4*1024;   % number of subfilters ~ FFT size
        else
          M0 = qLf*2*1024;   % number of subfilters ~ FFT size
        end

        for Mx = Mx_list
          i_case = i_case + 1;

          M = Mx*M0;

          fs_in = M0;
          fs = fs_in/M;
          Ts = 1/fs;
          n_out = t_max0/Ts;
          bin_bw = fs;
          N = n_out;
          n_in = (n_out+N_tap)*M;

          freq = [-Lf*M/2:Lf*M/2-1]/(Lf*M)*fs_in;
          n_freq = Lf*M;
          i_freq0 = Lf*M/2+1;
          ii_freq = i_freq0 + [floor(-2*Lf):ceil(2*Lf)];

          df_bin = freq(2)-freq(1);

          fprintf('fs=%.1f, Ntap=%.0f %s %.1fx\n',fs,N_tap,window_name,Lf);

          %
          % run PFB with noise only
          %

          sigma = sqrt(M0); A = 0; f1 = []; df_dt = [];
          [noise_out] = gen_chirp_pfb_matrix(...
                    sigma,A,f1,df_dt,fs,n_out,M,N_tap,window_name,Lf);
          bin_noise_db = 20*log10(rms(noise_out(:)))';

          snr_db = Inf;

          for i_chirp = 1:n_chirp
            m_chirp = m_chirp_all(i_chirp);
            f_offset = f_offset_all(i_chirp);
            df_dt = df_dt_all(i_chirp);

            fprintf('Case %.0f of %.0f, Chirp %.0f of %.0f\n',...
                  i_case,n_case,i_chirp,n_chirp);

            f1 = f_offset;
            f2 = f1 + df_dt*n_out*Ts;
            f_ctr = (f1+f2)/2;
            delta_f = f2 - f1;

            nb = max(8,ceil(m_chirp*Lf));
            nb2 = ceil(nb/2);

            config_str0 = sprintf('fs=%.1f Hz, T-max=%.0f',fs,t_max0);
            config_str1 = sprintf('Ntap=%.0f %s %.1fx, fs=%.1f',...
              N_tap,window_name,Lf,fs);
            %config_str1 = sprintf('Ntap=%.0f %s %.1fx, fs=%.1f, T-max=%.0f',...
            %  N_tap,window_name,Lf,fs,t_max0);

            config_str0_fn = sprintf('M-%.0f-fs-%.1f',...
              M,fs);  
            config_str1_fn = sprintf('%.2fx-Ntap-%.0f-%s-fs-%.1f',...
              Lf,N_tap,window_name,fs);

            %sig_str = sprintf('Chirp %.2f..%.2f Hz (%.3f Hz/sec)',...
            %  f1,f2,df_dt);
            sig_str = sprintf('%.2f Hz/sec',df_dt);
            sig_str_fn = sprintf('%.2fHzsec',df_dt);
            if ~isempty(mod_str)
              sig_str = sprintf('%s %s',sig_str,mod_str);
              sig_str_fn = [sig_str_fn '-' mod_str_fn];
            end

            sig_config_str0 = sprintf('%s %s',sig_str,config_str0);
            sig_config_str1 = sprintf('%s %s',sig_str,config_str1);

            sig_config_str0_fn = sprintf('%s-%s',sig_str_fn,config_str0_fn);
            sig_config_str1_fn = sprintf('%s-%s',sig_str_fn,config_str1_fn);

            %
            % generate input waveform and run PFB
            %

            fprintf('PFB %.1fx      ',Lf);

            tstart = tic;
            sigma = 0; A = 1;
            [x_out,freq1,t_out,coef] = gen_chirp_pfb_matrix(...
                      sigma,A,f1,df_dt,fs,n_out,M,N_tap,window_name,Lf,mod_params);
            fprintf(' %.1f sec\n',toc(tstart));

            t_max = n_in/fs_in;
            t_in = [0:n_in-1]/fs_in;

            xx = abs(x_out).^2;


            turbo_best_mean = mean(max(xx));

            turbo_best_db = 10*log10(turbo_best_mean);
            rel_snr_best_db = turbo_best_db - bin_noise_db;

            freq_t_out = interp1(t_in,f1+df_dt*t_in,t_out);

            x_out_db = 20*log10(max(1e-3,abs(x_out)));

            %
            % plot PFB bin output characteristics
            %

            if (i_chirp==1)
              plot_pfb_sg = 1;
            else
              plot_pfb_sg = 0;
            end
            plot_pfb_peaks = 1;

            plot_pfb_characteristics;

          end
        end

      end % window
    end  % N_tap
  end  % Lf
end  % modulation

