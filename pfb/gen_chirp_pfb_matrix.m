function [x_out,f_bin,t_out,coef,x_in] = gen_chirp_pfb_matrix(...
        sigma,A,f1,df_dt,fs,n_out,M,N_tap,window_name,Lf,mod_params,class_name,use_mex);
%
% Function to generate chirp waveform and run through polyphase or FFT filter bank 
%
% Notes
% Set sigma=0 for no noise.  sigma = sqrt(M) will provide approximately
%     unit amplitude noise at PFB output, but will depend on window type
% Set A, f1, df_dt to [] for no signal, or alternatively set A=0
% 
% inputs
%
% sigma        1x1        noise rms amplitude volts at PFB input {sqrt(I^2 + Q^2)}
% A            1xn_sig    chirp rms amplitude volts at PFB input {sqrt(I^2 + Q^2)}
% f1           1xn_sig    start frequency at t=0
% df_dt        1xn_sig    drift rate Hz/sec
% fs           1x1        sample rate at PFB output
% n_out        1x1        number of time samples
% M            1x1        Decimation factor = number of filter channels at 1x
% N_tap        1x1        Taps per subfilter, N_tap=1 for FFT only
% window_name  string     window type
% Lf           1x1        =1 for standard critically sampled PFB ("1x")
%                         =2 for 50% overlapped frequency bins ("2x")
%                         =4 for 50% overlapped frequency bins ("4x")
% mod_params   struct     modulation parameters, [] for none
% class_name   string     output type, ='double' (default), 'single'
% use_mex      1x1        =1 to use mex fns (default), =0 use m fns
%
% outputs
%
% x_out        n_freq x n_out complex voltage time waveforms for each frequency bin
% f_bin        n_freq x 1     bin center freq values
% t_out        1 x n_out      time sec for time series in each pfb bin
% coef         M x N_tap      PFB coefficients
% x_in         n_in x 1       complex voltage input time waveform (if needed)
%
% where n_freq = M (1x),2M (2x) and 4M (4x)%
%

  if (~exist('N_tap','var')),  N_tap = 8; end
  if isempty(N_tap),           N_tap = 8; end

  if (~exist('window_name','var')),  window_name = 'Hann'; end
  if isempty(window_name),           window_name = 'Hann'; end

  if (~exist('Lf','var')),  Lf = 1; end
  if isempty(Lf),           Lf = 1; end

  if (~exist('mod_params','var')),  mod_params = []; end
  if isempty(mod_params),           mod_params = []; end

  if (~exist('class_name','var')),  class_name = 'double'; end
  if isempty(class_name),           class_name = 'double'; end

  if (~exist('use_mex','var')),  use_mex = 1; end
  if isempty(use_mex),           use_mex = 1; end

  if mod(M,2)~=0
    error(sprintf('Error in gen_chirp_pfb_matrix, M must be even, %.0f\n',M));
  end

  if Lf<1
    error(sprintf('Error in gen_chirp_pfb_matrix, Lf must be >1, %.0f\n',Lf));
  end

  A = A(:)';
  f1 = f1(:)';
  df_dt = df_dt(:)';
  
  n_sig1 = length(A);
  n_sig2 = length(f1);
  n_sig3 = length(df_dt);
  
  if (max(A)==0)
    n_sig = 0;
  elseif (n_sig1==n_sig2) & (n_sig1==n_sig3)
    n_sig = n_sig1;
  else
    n_sig = max([n_sig1,n_sig2,n_sig3]);
    if (n_sig1==1)
      A = A*ones(1,n_sig);
    end
    if (n_sig2==1)
      f1 = f1*ones(1,n_sig);
    end
    if (n_sig3==1)
      df_dt = df_dt*ones(1,n_sig);
    end
    if (length(A)~=length(f1)) | (length(A)~=length(df_dt))
      error(sprintf(...
        'Error in gen_chirp_pfb_matrix, A, f1 and df_dt must be same length, %.0f vs %.0f vs %.0f\n',...
        n_sig1,n_sig2,n_sig3));
    end
  end

  fs_in = fs*M;
  
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
  % generate input waveform
  %

  fs_in = M*fs;
  n_in = M*n_out;

  if (sigma>0)
    x_in = sigma/sqrt(2)*complex(randn(n_in,1,class_name),randn(n_in,1,class_name));
  else
    x_in = complex(zeros(n_in,1,class_name),zeros(n_in,1,class_name));
  end
    
  if (n_sig>0)
    if isempty(mod_params)
      for i = 1:n_sig
        x_in = x_in + A(i)*gen_chirp_wf(fs_in,n_in,f1(i),df_dt(i),0,class_name);
      end
    else
      for i = 1:n_sig
        x_in = x_in + A(i)*gen_chirp_wf(fs_in,n_in,f1(i),df_dt(i),0,class_name)...
               .*gen_simple_mod(fs_in,n_in,mod_params,class_name);
      end
    end
  end
  
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




