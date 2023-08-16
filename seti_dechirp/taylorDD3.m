function [det_DD,freq2,df_dt_list,m_list] = ...
                 taylorDD3(xx,freq1,T_line,Lf,Lr,Nt,N0,df_dt_min,df_dt_max);
%
% Function to run energy detection with input xx and 
% a range of Doppler drift rates based on a Taylor tree algorithm
% Updated algorithm that avoids bit-reversed indexing but produces identical
% results to legacy taylor_flt C Code
% 
% Outputs
%
% inputs
%
% xx           n_freq1 x n_time  mag squared time waveforms for each frequency bin (power)
% freq1        n_freq1 x 1       PFB bin center freq values
% T_line       1 x 1             time sec for each line in spectrogram
% Lf           1 x 1             overlap factor in spectrogram (power of 2)
% Lr           1 x 1             upsample for drift rate increment
% Nt            1 x 1             last stage number of time samples (Nt/N0 is power of 2)
% N0           1 x 1             first stage number of time samples
% df_dt_min    1 x 1             minimum drift rate to be computed >= -fs^2 (approx)
% df_dt_max    1 x 1             maximum drift rate to be computed <= fs^2 (approx)
%
% outputs
%
% det_DD       n_freq2 x n_m     integrated energy array
% freq2        n_freq2 x 1       det_DD bin center freq values
% df_dt_list   n_m x 1           det_DD drift rate Hz/sec array
% m_list       n_m x 1           drift rate index array
%
% where
%
% n_freq2 = n_freq1/max(Lf/Lr,1)
% n_m = depends on df_dt_min, df_dt_max
%
% Notes
%
% class of det_DD will be same as class of input xx, e.g. 'double' or 'single'
%
% Normally n_time will equal Nt.  If n_time<Nt, input will be zero padded and
% only n_time lines will be summed. If n_time>Nt input will be truncated and
% only Nt integrations will be performed.
%
% Define bin_bw = bin bandwidth in each pfb bin = fs = sampling rate at
% bin output
%
% Input in xx will have frequency oversampling factor Lf, so 
%    delta_freq1 = bin_bw/Lf = freq1(2) - freq1(1)
% Output in det_DD does critical frequency sampling factor, so 
%    delta_freq2 = bin_bw/Lr = freq2(2) - freq2(1)
% Note Lf >= Lr, and Lf must be an integer multiple of Lr
%
% Output will also have a df_dt increment
%    delta_df_dt = df_dt_list(2) - df_dt_list(1),
% where delta_df_dt = fs/(Lr*n_time*T_line)
%
% Normally T_line will equal 1/fs.  Sometimes the spectrogram lines will be
% averaged N_preDD times prior to input, so that T_line = N_preDD/fs.
% However, this is generally not recommended (though accommodated here),
% because the frequency drifts over fs^2/N_preDD will suffer significant
% attenuation in the averaging process.
% 

  if (~exist('N0','var')),  N0 = 4; end
  if isempty(N0),           N0 = 4; end

  [n_freq1,n_time] = size(xx);
  if length(freq1)~=n_freq1
    error(sprintf('Error in taylorDD3, xx and freq1 sizes incompatible, %.0f vs %.0f\n',...
      n_freq1,length(freq1)));
  end

  if Nt~=N0*2.^round(log2(Nt/N0))
    error(sprintf('Error in taylorDD3, Nt/N0 not a power of 2, %.0f\n',Nt));
  end
  
  if (Lr>Lf)
    fprintf('Warning in taylorDD3, Lr=%.0f cannot be greater than Lf=%.0f\n',...
      Lr,Lf);
    Lr = Lf;
  end
  if mod(Lf,Lr)~=0
    error(sprintf('Error in taylorDD3, Lf=%.0f must be multiple of Lr=%.0f\n',...
      Lf,Lr));
  end
  n_freq2 = n_freq1*Lr/Lf;

  df_bin = freq1(2)-freq1(1);  % may be less than 1/Ts if bins are overlapped
  fs = df_bin*Lf;
  freq2 = freq1(1:(Lf/Lr):end);
  
  if (~exist('df_dt_min','var')),  df_dt_min = 0; end
  if isempty(df_dt_min),           df_dt_min = 0; end

  if (~exist('df_dt_max','var')),  df_dt_max = fs.^2; end
  if isempty(df_dt_max),           df_dt_max = fs.^2; end

  n_stage = log2(Nt/N0)+1;
  
  m_max0 = ceil(df_dt_max*T_line/fs*(Lr*N0-1));
  m_min0 = floor(df_dt_min*T_line/fs*(Lr*N0-1));
  [m_min,m_max] = m_limit(m_min0,m_max0,Nt,N0);  % limits at output
  m_list  = [m_min:m_max];
  delta_df_dt = fs/T_line/(Nt*Lr-1);
  df_dt_list = m_list*delta_df_dt;

  dd0 = max(Nt*Lf,ceil(Nt*Lf*max(abs(df_dt_list))*T_line/fs));
  dd0 = 2.^ceil(log2(dd0));

  Nr0 = m_max0 - m_min0 + 1;
  Nr  = m_max  - m_min  + 1;

  %
  % set up energy work array
  %
  
  if (n_time<Nt)
    xx = [xx zeros(n_freq1,Nt-n_time)];
  end
  if (n_freq1<Nt)
    xx = [xx ; zeros(Nt - n_freq1,Nt)];
  end

  % first stage processing

  det_DD_i = taylor_stage1(xx,n_freq1,Lf,Lr,Nt,N0,m_min0,m_max0,dd0);
  
  % i-th stage processing

  % extended Taylor fast stages
  for i_stage = 2:n_stage
    det_DD_i = taylor_stage_i(det_DD_i,i_stage,n_freq2,Nt,N0,m_min0,m_max0,dd0);
  end

  %
  % final result
  %
  
  det_DD = det_DD_i(dd0+1:dd0+n_freq2,:)/n_time;
  
end

function det_DD1 = taylor_stage1(xx,n_freq1,Lf,Lr,Nt,N0,m_min0,m_max0,dd0)
  %
  % "Slow" first stage
  %
  Nr0 = m_max0 - m_min0 + 1;
  freq_bin_offset = Lf/Lr*[m_min0:m_max0]'*[0:N0-1]/(N0-1);
  int_bin_offset = round(freq_bin_offset);
  
  % augment xx with zeros

  class_name = class(xx);
  xx = [zeros(dd0,Nt,class_name) ; xx ; zeros(dd0,Nt,class_name)];
    
  n_group = Nt/N0;  % number of butterflies

  n_freq2 = n_freq1*Lr/Lf;
  det_DD1 = zeros(n_freq2+2*dd0,n_group*Nr0,class_name);  

  for ig = 1:n_group
    it0 = (ig-1)*N0;
    m0  = (ig-1)*Nr0;
    b_sum = zeros(n_freq2,Nr0,class_name);
    for mm = 1:Nr0
      for it = 1:N0
        k_ofs = dd0 + int_bin_offset(mm,it);
        b_sum(:,mm) = b_sum(:,mm) + xx(k_ofs+1:(Lf/Lr):k_ofs+n_freq1,it0+it);
      end
    end
    det_DD1(dd0+1:dd0+n_freq2,m0+[1:Nr0]) = b_sum;
  end
end

function det_DD_out = taylor_stage_i(det_DD_in,i_stage,n_freq2,Nt,N0,m_min0,m_max0,dd0)

  N_out = N0* 2.^(i_stage-1);  % number of time samples integrated in this stage
  N_in = N_out/2;                % number of time samples integrated in previous stage
  n_group = Nt/N_out;
  
  [m_min_out,m_max_out] = m_limit(m_min0,m_max0,N_out,N0); % drift rate index limits - this stage
  Nr_out = m_max_out - m_min_out + 1;
  [m_min_in,m_max_in] = m_limit(m_min0,m_max0,N_in,N0); % drift rate index limits - previous stage
  Nr_in = m_max_in - m_min_in + 1;
  
  class_name = class(det_DD_in);
  det_DD_out = zeros(n_freq2+2*dd0,n_group*Nr_out,class_name);
  
  for ig = 1:n_group
    m00 = (ig-1)*2*Nr_in;
    m01 = m00 + Nr_in;
    m02 = (ig-1)*Nr_out;

    for m = m_min_out:m_max_out
      i2 = m - m_min_out;
      if (m>=0)
        k_ofs = dd0 + floor((m+1)/2);
        mby2 = floor(m/2);
      else
        k_ofs = dd0 - floor((-m+1)/2);
        mby2 = -floor(-m/2);
      end
      i1 = mby2 - m_min_in;
      det_DD_out(dd0+1:dd0+n_freq2,m02+i2+1) = ...
                             det_DD_in(dd0+1:dd0+n_freq2,m00+i1+1) + ...
                             det_DD_in(k_ofs+1:k_ofs+n_freq2,m01+i1+1);
    end    
  end
end

function [m_min_i,m_max_i] = m_limit(m_min0,m_max0,Ni,N0)
  % index limits for Taylor DD
  if (m_max0>0)
    m_max_i = (m_max0+1)*Ni/N0 - 1;
  else
    m_max_i = m_max0*Ni/N0;
  end
  if (m_min0>=0)
    m_min_i =  m_min0*Ni/N0;
  else
    m_min_i = (m_min0-1)*Ni/N0 + 1;
  end
  m_max_i = round(m_max_i);
  m_min_i = round(m_min_i);
end

