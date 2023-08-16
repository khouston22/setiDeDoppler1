%
% script to verify legacy TurboSETI gives identical results to 
% matlab-equivalent function turbo_seti1.m (no bitrev() use)
% also can check equivalence to turbo_seti0_py.m and turbo_seti0_c.m
% which run legacy code within matlab functions
%
% Note: to run Python functions successfully using Anaconda Python install, 
% may need to start Matlab from Anaconda Prompt window, e.g. in Windows
% cd C:\Program Files\MATLAB\R2017b\bin
% C:\Program Files\MATLAB\R2017b\bin>matlab
%

clear;
%clear classes;

addpath('..\seti_dechirp');

test_python = 0;    % =1  use python version, =0 use taylor_flt.c

use_random = 1;     % =1  use random input data, =0 create drift lines
pos_neg = 0;        % =1 positive drift rates only, =0 pos and negative
print_matrices = 0; % =1 print mex and .m outputs and differences

if (0)
  Nt_list = [16 32 64 128 256 512 1024];
elseif (0)
  Nt_list = [16];
else 
  Nt_list = [1024];
end

Lf_list = [1 3/2 2];

if test_python
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
  
  if (0)
    for nbits=1:7
      ibitrev = [];
      for i=[0:2.^nbits - 1]
        ibitrev(i+1) = py.bitrev1.bitrev(uint64(i),uint64(nbits));
      end
      if  (2.^nbits > 16), ibitrev=reshape(ibitrev,16,[])'; end
      print_compact_matrix(sprintf('Python bitrev test, nbits=%.0f',nbits),ibitrev,3,0);
    end
  end
else
  mex taylor_flt.c
end

%return;

for Nt = Nt_list
  
   n_freq=8*Nt; 
  
  if (use_random)
    m0_list = [0 0 0];
  elseif (0)
    if (pos_neg==1)
      m0_list = [0 Nt/8 Nt/4 Nt/2 Nt];
    else
      m0_list = [0 Nt/8 Nt/4 Nt/2 Nt];
      m0_list = [-m0_list(end:-1:2) m0_list]; 
    end
  else
    m0_list = [Nt/8];
    %m0_list = [-Nt/8];
  end

  for Lf = Lf_list
    Lr = Lf;
    fs = 1;
    df_bin = fs/Lf;
    n_freq=8*Nt*Lf; 

    for m0 = m0_list
      in = zeros(n_freq,Nt,'single');
      if (use_random)
        %in = round(randn(n_freq,Nt));
        in = abs(round(10*randn(n_freq,Nt)));
        k0=NaN;
      else
        k0=Nt+8;
        for i=1:Nt
          if (m0==0)
            in(k0+1,i)=11;                 % m0 = 0
          else
            denom = Nt/m0;
            in(k0+1+round((i-1)/denom),i)=11;  % m0 = Nt/denom
          end
        end
      end

      % legacy taylor_flt
      n_freq1=n_freq+2*Nt; 
      mlen = n_freq1*Nt; 
      in1 = [in ; zeros(2*Nt,Nt)];
      inbuf = single(reshape(in1,mlen,[]));

      t0 = tic;
      if test_python
        py_inbuf1 = py.numpy.array(inbuf(:).');
        py_outbuf1 = py.taylor_flt.flt(py_inbuf1,uint64(mlen),uint64(Nt));
        out1a = ndarray2dbl(py_outbuf1,n_freq1)/Nt;
        out1 = bitrev_cols(out1a);
      else
        outbuf1 = taylor_flt(inbuf,mlen,Nt)';
        out1 = reshape(outbuf1,n_freq1,[])/Nt;
      end

      if (pos_neg==1)
        out = out1(1:n_freq,:);
      else
        in2 = [zeros(2*Nt,Nt); in];
        inbuf = single(reshape(flipud(in2),mlen,[]));
        if test_python
          py_inbuf2 = py.numpy.array(inbuf(:).');
          py_outbuf2 = py.taylor_flt.flt(py_inbuf2,uint64(mlen),uint64(Nt));
          out2a = ndarray2dbl(py_outbuf2,n_freq1)/Nt;
          out2 = bitrev_cols(out2a);
        else
          outbuf2 = taylor_flt(inbuf,mlen,Nt)';
          out2 = reshape(outbuf2,n_freq1,[])/Nt;
        end
        out2 = fliplr(flipud(out2));

        out = [out2(2*Nt+1:end,1:end-1) out1(1:n_freq,:)];
      end
      t_elapsed = toc(t0);

      fs = 1;
      df_bin = fs;
      f_pfb = [-n_freq/2:n_freq/2-1]*df_bin;
      t_pfb = [0:Nt-1]'/fs;
      N0=4;
      Lf = 1;
      Lr = 1;
      T_line = 1/fs;
      if (0)
        % new taylor alg without bitrev()
        [det_DD,freq2,df_dt_list,m_list] = ...
                            taylorDD1(in,f_pfb,t_pfb,Nt,N0,pos_neg);
      elseif (0)
        % upgraded taylor alg with Lf Lr   
        df_dt_min = -fs.^2;
        df_dt_max =  fs.^2;
        [det_DD,freq2,df_dt_list,m_list] = ...
                     taylorDD3(in,f_pfb,T_line,Lf,Lr,Nt,N0,df_dt_min,df_dt_max);
      elseif (1)
        % upgraded taylor alg within fastDD  
        df_dt_min = -fs.^2;
        df_dt_max =  fs.^2;
        [det_DD,freq2,df_dt_list,m_list] = ...
                     fastDD3(in,f_pfb,T_line,Lf,Lr,Nt,N0,df_dt_min,df_dt_max,0);
      else
        %
        % check equivalent m functions for consistency
        %
        if test_python
          % function using _taylor_tree/_core_numba.py
          [det_DD,df_dt_list,m_list] = taylorDD0_py(in,f_pfb,t_pfb,Nt,pos_neg);
        else
          % function using legacy taylor_flt.c
          [det_DD,df_dt_list,m_list] = taylorDD0_c(in,f_pfb,t_pfb,Nt,pos_neg);
        end
      end

      Nr = length(m_list);
      difference = out - det_DD;

      if (print_matrices)
        print_compact_matrix('in',in',2,0);
        %print_compact_matrix('flipud(in)',flipud(in)',2,0);
        print_compact_matrix('out',out',2,0);
        %print_compact_matrix('out1',out1',2,0);
        %print_compact_matrix('out2',out2',2,0);
        print_compact_matrix('det_DD',det_DD',2,0);
        print_compact_matrix('difference',difference',2,0);
      end

      fprintf('Nt=%4.0f Nr=%.0f Nf=%4.0f k0=%4.0f m0=%4.0f max diff=%2.0f, %5.3f sec\n\n',...
        Nt,Nr,n_freq,k0,m0,max(max(abs(difference))),t_elapsed);
    end
  end
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


