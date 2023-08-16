%
% script to verify C-code mex version of pfb_Lfx_mex.c gives identical results to 
% matlab function gen_chirp_pfb_matrix.m and gen_chirp_pfb_matrix2.m 
%

clear;

addpath('..\misc_fns');
addpath('..\pfb');

use_random = 1;     % =1  use random input data, =0 create drift lines
print_matrices = 0; % =1 print mex and .m outputs and differences
do_fft = 1;         % =1 do c fft calls (vs. ffts done in matlab)

if (1)
  % delete pfb_Lfx_mex.mexw64
  fprintf('Testing pfb_Lfx_mex.c, do_fft=%.0f\n',do_fft);
  mex pfb_Lfx_mex.c pfb_fns.c kiss_fft.c
end

%return;

if (1)
  %N_tap_list = [1,3,4,6,8];
  N_tap_list = [4];
  n_trial = 1;
  Lf_list = [1 1 5/4 4/3 3/2 7/4 2 3 4];
  n_out = 512;
%   window1_list = {'RectWin','Hann'};
%   window2_list = {'Hann','lpfxxx'};
  window1_list = {'RectWin'};
  window2_list = {'lpfxxx'};
elseif (0)
  N_tap_list = [4];
  n_trial = 4;
  Lf_list = [1 3/2 2 4];
  n_out = 512;
  window1_list = {'RectWin'};
  window2_list = {'lpfxxx'};
else
  N_tap_list = [3];
  n_trial = 1;
  Lf_list = [3/2];
  n_out = 8;
  window1_list = {'RectWin'};
  window2_list = {'Hann'};
  print_matrices = 1;
end


for Lf = Lf_list
  for i_trial = 1:n_trial
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
          M = 8*n_out;   % number of subfilters = decimation
        else
          M = 8*round(n_out/qLf)*qLf;   % number of subfilters ~ FFT size
        end
        fs_in = 1.0*M;

        n_in = (n_out+N_tap)*M;
        %n_in = n_out*M;

        fs = fs_in/M;
        t_max = n_out/fs;

        freq = [-Lf*M/2:Lf*M/2-1]'/(Lf*M)*fs_in;
        n_freq = Lf*M;

        if (N_tap>2)
          apply_sinc = 1;
        else
          apply_sinc = 0;
        end

        coef = single(calc_sinc_window_coefs(M,N_tap,window_name,apply_sinc));

        %
        % input signal = noise
        %

        x_in = complex(round(10*randn(n_in,1,'single')),round(10*randn(n_in,1,'single')));

        %
        % run pfb or FFT filter bank - m file
        %

        t_start = tic;
        [x_out_m,fs_out,f_bin,t_out,xcoef] = run_pfb_Lfx(x_in,coef,fs_in,n_out,pLf,qLf);
        t_elapsed_m = toc(t_start);

        %
        % run pfb or FFT filter bank - mex file
        %
        t_start = tic;
        [x_out_mex,fs_out_mex,f_bin_mex,t_out_mex,xcoef_mex] = ...
                           pfb_Lfx_mex(x_in,coef,fs_in,n_out,pLf,qLf,do_fft);
        if (do_fft==0)
          % do FFTs and interleave
          N_fft = n_freq/pLf;
          for i_out = 1:n_out
            temp = x_out_mex(:,i_out);
            for i1=1:pLf
              temp2 = temp((i1-1)*N_fft + [1:N_fft]);
              temp2 = fftshift(fft(temp2,N_fft));
              x_out_mex(i1:pLf:n_freq,i_out) = temp2;
            end
          end
        end
        t_elapsed_mex = toc(t_start);

        difference = x_out_mex - x_out_m;
        max_diff = max(max(abs(difference)));

        if (print_matrices)
          print_compact_matrix('real x_in',real(x_in.'),3,0);
          print_compact_matrix('imag x_in',imag(x_in.'),3,0);
          print_compact_matrix('real x_out_mex',real(x_out_mex.'),3,0);
          print_compact_matrix('imag x_out_mex',imag(x_out_mex.'),3,0);
          print_compact_matrix('real x_out_m',real(x_out_m.'),3,0);
          print_compact_matrix('imag x_out_m',imag(x_out_m.'),3,0);
          print_compact_matrix('difference',difference',3,0);
        end

%         if (1 && pLf>1)
%           print_compact_matrix('real xcoef*M 1',real(M*xcoef(:,:,1).'),3,1);
%           print_compact_matrix('real xcoef_mex 1',real(M*xcoef_mex(:,:,1).'),3,1);
%           print_compact_matrix('difference 1',M*xcoef_mex(:,:,1).'-M*xcoef(:,:,1).',3,1);
%           print_compact_matrix('real xcoef*M 2',real(M*xcoef(:,:,2).'),3,1);
%           print_compact_matrix('real xcoef_mex 2',real(M*xcoef_mex(:,:,2).'),3,1);
%           print_compact_matrix('difference 2',M*xcoef_mex(:,:,2).'-M*xcoef(:,:,2).',3,1);
%           print_compact_matrix('xcoef 2 degrees',180/pi*angle(xcoef(:,:,2).'),4,0);
%           print_compact_matrix('xcoef_mex 2 degrees',180/pi*angle(xcoef_mex(:,:,2).'),4,0);
%           print_compact_matrix('difference 2 degrees',180/pi*(angle(xcoef_mex(:,:,2).')-angle(xcoef(:,:,2).')),4,0);
%         end
        
        fprintf('Lf=%4.2f N_tap=%2.0f M=%4.0f %7s n_freq=%5.0f n_out=%3.0f max diff=%2.0f, %4.0f mex vs %4.0f m msec\n\n',...
        Lf,N_tap,M,window_name,n_freq,n_out,max_diff,t_elapsed_mex*1e3,t_elapsed_m*1e3);

      end
    end
  end
end


