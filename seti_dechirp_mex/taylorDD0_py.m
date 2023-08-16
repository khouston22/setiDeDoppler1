function [det_DD,df_dt_list,m_list] = taylorDD0_py(xx,freq1,t_pfb,N,pos_neg);
%
% Function to run energy detection with input xx and 
% a range of Doppler drift rates based on legacy taylor_flt C Code
% Note: can operate on 1x frequency sampling (Lf=1) only
% 
% Outputs
%
% inputs
%
% xx           n_freq1 x n_time  max squared time waveforms for each frequency bin (power)
% freq1        n_freq1 x 1       PFB bin center freq values
% t_pfb        n_time x 1        time sec for time series in each pfb bin
% N            1 x 1             last  stage number of time samples (power of 2)
% pos_neg      1 x 1             =1 for positive drift rates (m=0..n_time-1)
%                                =-1 for negative drift rates (m=-n_time+1..0)
%                                =0 for pos & neg drift rates (m=-n_time+1..n_time-1) (default)
%
% outputs
%
% det_DD       n_freq1 x n_m     integrated energy array
% df_dt_list   n_m x 1           drift rate Hz/sec array
% m_list       n_m x 1           drift rate index array
%

  if (~exist('pos_neg','var')),  pos_neg = 0; end
  if isempty(pos_neg),           pos_neg = 0; end

  [n_freq1,n_time] = size(xx);
  %N = 2.^round(log2(n_time));
  
  if length(freq1)~=n_freq1
    error(sprintf('Error in taylorDD0, xx and freq1 sizes incompatible, %.0f vs %.0f\n',...
      n_freq1,length(freq1)));
  end
  if (length(t_pfb)~=n_time) || (N~=n_time)
    error(sprintf('Error in taylorDD0, xx, t_pfb sizes incompatible with N=%.0f, %.0f vs %.0f\n',...
      n_time,length(t_pfb)),N);
  end

  if N~=2.^round(log2(N))
    error(sprintf('Error in taylorDD0, N not a power of 2, %.0f\n',N));
  end
  
  [v, e, loaded] = pyversion;
  if ~loaded
    pyversion C:\Users\KMH\anaconda3\python.exe
  end
%   [own_path, ~, ~] = fileparts(mfilename('fullpath'));
%   module_path = fullfile(own_path, '..');
%   python_path = py.sys.path;
%   if count(python_path, module_path) == 0
%       insert(python_path, int32(0), module_path);
%   end
  py.importlib.import_module('numpy');
  py.importlib.import_module('numba');
  py.importlib.import_module('bitrev1');
  py.importlib.import_module('taylor_flt');

  
  df_bin = freq1(2)-freq1(1);  % may be less than 1/Ts if bins are overlapped
  Ts = t_pfb(2)-t_pfb(1);
  Ts_df_bin = Ts*df_bin;
  
  if (abs(Ts_df_bin-1)<.001)
    Lf = 1;
  else
    error(sprintf('Error in taylorDD0, df_bin=%.3f x Ts=%.3f=%.3f must be 1\n',...
      df_bin,Ts,Ts_df_bin));
  end

  if (pos_neg==1)
    m_list  = [0:N-1];
  elseif (pos_neg==-1)
    m_list  = [-(N-1):0];
  else
    m_list  = [-(N-1):N-1];
  end

  fs = 1/Ts;
  df_dt_list = m_list*fs.^2/(N-1);
  
  %
  % set up energy work array
  %

  if (n_time<N)
    xx = [xx zeros(n_freq1,N-n_time)];
  end
  if (n_freq1<N)
    xx = [xx ; zeros(N - n_freq1,N)];
  end
    
  %
  % processing stages
  %
  
  out_p = single([]);
  out_n = single([]);

  if (pos_neg>=0)
    %
    % positive drift rates
    % all processing stages
    %

    n_freq1x=n_freq1+2*N; 
    mlen = n_freq1x*N; 
    inbuf = single(reshape([xx ; zeros(2*N,N)],mlen,[]));
    
    py_inbuf1 = py.numpy.array(inbuf(:).');
    py_outbuf1 = py.taylor_flt.flt(py_inbuf1,uint64(mlen),uint64(N));
    out_p = ndarray2dbl(py_outbuf1,n_freq1x)/N;
    out_p = bitrev_cols(out_p);
    out_p = out_p(1:n_freq1,:);
  end

  if (pos_neg<=0)
    %
    % negative drift rates
    % all processing stages
    %
    
    % all processing stages

    n_freq1x=n_freq1+2*N; 
    mlen = n_freq1x*N; 
    inbuf = single(reshape(flipud([zeros(2*N,N) ; xx]),mlen,[]));
    
    py_inbuf2 = py.numpy.array(inbuf(:).');
    py_outbuf2 = py.taylor_flt.flt(py_inbuf2,uint64(mlen),uint64(N));
    out_n = ndarray2dbl(py_outbuf2,n_freq1x)/N;
    out_n = bitrev_cols(out_n);
    out_n = fliplr(flipud(out_n));
    out_n = out_n(2*N+1:end,:);

    if (pos_neg==0)
      out_n = out_n(:,1:end-1);
    end
  end

  %
  % final result
  %

  det_DD = double([out_n out_p]);
  
end

function out = ndarray2dbl(py_outbuf,nx)
  pyList = py_outbuf.tolist();
  out = cellfun(@double,cell(pyList));
  out = reshape(out,nx,[]);
end

function out =  bitrev_cols(in)
  [m,n] = size(in);
  out = NaN(m,n);
  nbits = log2(n);
  for i=[0:n-1]
    ibitrev1 = py.bitrev1.bitrev(uint64(i),uint64(nbits));
    out(:,ibitrev1+1) = in(:,i+1);
  end
end

