function [x_out,f_bin,t_out,coef] = gen_chirp_pfb_matrix2(...
               x_in,fs,n_out,M,N_tap,window_name,Lf,use_mex);
%
% Function to run a waveform and run through polyphase or FFT filter bank 
% 
% inputs
%
% x_in         n_in x 1   complex voltage input time waveform (if needed)
% fs           1x1        sample rate at PFB output
% n_out        1x1        number of time samples
% M            1x1        Decimation factor = number of filter channels at 1x
% N_tap        1x1        Taps per subfilter, N_tap=1 for FFT only
% window_name  string     window type
% Lf           1x1        =1 for standard critically sampled PFB ("1x")
%                         =2 for 50% overlapped frequency bins ("2x")
%                         =4 for 50% overlapped frequency bins ("4x")
% use_mex      1x1        =1 to use mex fns (default), =0 use m fns
%
% outputs
%
% x_out        n_freq x n_out complex voltage time waveforms for each frequency bin
% f_bin        n_freq x 1     bin center freq values
% t_out        1 x n_out      time sec for time series in each pfb bin
% coef         M x N_tap      PFB coefficients
%
% where n_freq = M (1x),2M (2x) and 4M (4x)%
%

  if (~exist('N_tap','var')),  N_tap = 8; end
  if isempty(N_tap),           N_tap = 8; end

  if (~exist('window_name','var')),  window_name = 'Hann'; end
  if isempty(window_name),           window_name = 'Hann'; end

  if (~exist('Lf','var')),  Lf = 1; end
  if isempty(Lf),           Lf = 1; end

  if (~exist('use_mex','var')),  use_mex = 1; end
  if isempty(use_mex),           use_mex = 1; end

  if mod(M,2)~=0
    error(sprintf('Error in gen_chirp_pfb_matrix2, M must be even, %.0f\n',M));
  end

  if Lf<1
    error(sprintf('Error in gen_chirp_pfb_matrix, Lf must be >1, %.0f\n',Lf));
  end

  fs_in = fs*M;
  n_in = M*n_out;
    
  n_freq = Lf*M;
  f_bin = [-Lf*M/2:Lf*M/2-1]/(Lf*M)*fs_in;

  %
  % generate pfb or window coefficients
  %
  
  if (N_tap>2)
    apply_sinc = 1;
  else
    apply_sinc = 0;
  end
  
  coef = calc_sinc_window_coefs(M,N_tap,window_name,apply_sinc);
    
  %
  % run pfb or FFT filter bank
  %

  [pLf,qLf]=rat(Lf);

  if use_mex
    do_fft = 1;
    [x_out,fs_out,f_bin,t_out] = pfb_Lfx_mex(single(x_in),single(coef),fs_in,n_out,pLf,qLf,do_fft);
  else
    [x_out,fs_out,f_bin,t_out] = run_pfb_Lfx(x_in,coef,fs_in,n_out,pLf,qLf);
  end
end




