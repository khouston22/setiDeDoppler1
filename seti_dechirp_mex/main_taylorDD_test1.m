%
% script to verify C-code mex version of taylor_DD.c gives identical results to 
% matlab function taylorDD3.m 
% Also checks fast_DD.c against fastDD3.m
%

clear;
%clear classes;

addpath('..\seti_dechirp');
addpath('..\misc_fns');
addpath('..\pfb');

test_taylor = 1;     % =1  test taylorDD, =0 test fastDD

for test_taylor = [1 0]
  if test_taylor
    delete taylor_DD_mex.mexw64
    fprintf('Testing taylor_DD_mex.c\n');
    mex taylor_DD_mex.c taylor_DD_fns.c
  else
    delete fast_DD_mex.mexw64
    fprintf('Testing fast_DD_mex.c\n');
    mex fast_DD_mex.c taylor_DD_fns.c fast_DD_fns.c
  end

  %return;

  use_random = 1;     % =1  use random input data, =0 create drift lines
  pos_neg = 0;        % =1 positive drift rates only, =0 pos and negative
  print_matrices = 0; % =1 print mex and .m outputs and differences

  n_trial = 1;
  Lf_list = [1 3/2 2];
  Lr_list = Lf_list;
  alg_ID = 3;         % fastDD algorithm 0-6
  fprintf('alg_ID = %.0f\n',alg_ID);

  if (1)
    Nt_list = [16 32 64 128 256 512 1024];
    Lf_list = [1];
    n_trial = 4;
  elseif (1)
    Nt_list = [16 32 64 128 256 512 1024];
    n_trial = 1;
    Lf_list = [2];
    Lr_list = [1];
  else 
    Nt_list = [1024];
    n_trial = 8;
    Lf_list = [2];
  end

  for Nt = Nt_list

    N0=Nt;  % single stage only
    %N0=Nt/2; % 2 stages
    %N0=8;
    N0=4;

    if (use_random)
      m0_list = zeros(1,n_trial);
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

    for iLf = 1:length(Lf_list) 
      Lf = Lf_list(iLf);
      if length(Lr_list)==length(Lf_list)
        Lr = Lr_list(iLf);
      else
        Lr = Lf;
      end
      fs = 1;
      df_bin = fs/Lf;
      n_freq=8*Nt*Lf; 

      for m0 = m0_list
        in = zeros(n_freq,Nt,'single');
        if (use_random)
          in = abs(round(10*randn(n_freq,Nt,'single')));
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

        f_pfb = [-n_freq/2:n_freq/2-1]'*df_bin;
        t_pfb = [0:Nt-1]/fs;
        T_line = 1/fs;
        if (pos_neg==0)
          df_dt_min = -fs.^2;
          df_dt_max =  fs.^2;
        elseif (pos_neg==1)
          df_dt_min = 0;
          df_dt_max =  fs.^2;
        elseif (pos_neg==-1)
          df_dt_min = -fs.^2;
          df_dt_max =  0;
        end

        if test_taylor
          % upgraded taylor alg taylorDD3.m   
          t0 = tic;
          [det_DD,freq2,df_dt_list,m_list] = ...
                         taylorDD3(in,f_pfb,T_line,Lf,Lr,Nt,N0,df_dt_min,df_dt_max);
          t_elapsed_m = toc(t0);

          % mex taylor_DD
          t0 = tic;
          [det_DD_mex,freq2_mex,df_dt_list_mex,m_list_mex] = ...
                       taylor_DD_mex(in,f_pfb,T_line,Lf,Lr,Nt,N0,df_dt_min,df_dt_max);
          t_elapsed_mex = toc(t0);
        else
          % fastDD alg fastDD3.m   
          t0 = tic;
          [det_DD,freq2,df_dt_list,m_list] = ...
                         fastDD3(in,f_pfb,T_line,Lf,Lr,Nt,N0,df_dt_min,df_dt_max,alg_ID);
          t_elapsed_m = toc(t0);

          % mex fastDD
          t0 = tic;
          [det_DD_mex,freq2_mex,df_dt_list_mex,m_list_mex] = ...
                       fast_DD_mex(in,f_pfb,T_line,Lf,Lr,Nt,N0,df_dt_min,df_dt_max,alg_ID);
          t_elapsed_mex = toc(t0);
        end

         Nr = length(m_list);
         difference = det_DD_mex - det_DD;
         max_diff = max(max(abs(difference)));

        if (print_matrices)
          print_compact_matrix('in',in',2,0);
          print_compact_matrix('det_DD_mex',det_DD_mex',2,0);
          print_compact_matrix('det_DD',det_DD',2,0);
          print_compact_matrix('difference',difference',2,0);
        end

        fprintf('Nt=%4.0f Nr=%4.0f Nf=%4.0f max diff=%2.0f, %5.0f mex vs %5.0f m msec\n',...
          Nt,Nr,n_freq,max_diff,t_elapsed_mex*1e3,t_elapsed_m*1e3);
      end
    end
  end
end


