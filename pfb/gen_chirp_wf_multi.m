function x_in = gen_chirp_wf_multi(sigma,A,f1,df_dt,fs,n_out,M,mod_params,class_name);
%
% Function to generate chirp waveform that will subsequently be run through 
% polyphase or FFT filter bank 
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
% mod_params   struct     modulation parameters, [] for none
% class_name   string     output type, ='double' (default), 'single'
%
% outputs
%
% x_in         n_in x 1       complex voltage input time waveform (if needed)
%
%

  if (~exist('mod_params','var')),  mod_params = []; end
  if isempty(mod_params),           mod_params = []; end

  if (~exist('class_name','var')),  class_name = 'double'; end
  if isempty(class_name),           class_name = 'double'; end

  if mod(M,2)~=0
    error(sprintf('Error in gen_chirp_pfb_matrix, M must be even, %.0f\n',M));
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
  
end




