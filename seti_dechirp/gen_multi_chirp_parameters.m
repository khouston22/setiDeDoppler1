function [f1_all,f2_all,df_dt_all,m0_chirp_all,f_offset_all,f_total] = ...
          gen_multi_chirp_parameters(...
               fs,Nt_pfb,M,f_incr_offset,m0_chirp_list,m0_offset_list,f_offset_list);
%
% Function to set parameters for a multiplexed chirp waveform 
% that will be generated by gen_chirp_wf_multi() and then run through
% polyphase or FFT filter bank 
% 
% inputs
%
% fs             1x1        sample rate at PFB output
% Nt_pfb         1x1        number of time samples at PFB output
% M              1x1        Decimation factor = number of filter channels at 1x
% f_incr_offset  1x1        frequency increment between successive chirps
% m0_chirp_list  1xn_sig    nominal chirp drift rate index list
% m0_offset_list 1xn_m_ofs  list of offsets to drift rate index
% f_offset_list  1xn_f_ofs  list of frequency offsets
% 
% Note: if m0_offset_list==NaN, random value applied -.5 to +.5 
%       if f_offset_list==NaN, random value applied (-.5 to +.5)*bin_bw
%
% outputs
%
% f1_all        1xn_all    start frequency array at t=0 Hz
% f2_all        1xn_all    end frequency array at t=Nt_pfb/fs Hz
% df_dt_all     1xn_all    drift rate array Hz/sec
% m0_chirp_all   1xn_all    drift rate array index
% f_offset_all  1xn_all    initial offset array Hz
%
% where n_all = n_sig * n_m_ofs * n_f_ofs
%

  fs_in = fs*M;
  n_in = M*Nt_pfb;
    
  freq = [-M/2:M/2-1]/M*fs_in;
  t_max = Nt_pfb/fs;
  bin_bw = fs_in/M;
  
  n_chirp = length(m0_chirp_list)*length(m0_offset_list)*length(f_offset_list);
  
  %
  % determine range of chirp values
  %
  
  if any(isnan(f_offset_list))
    f_chirp_offset = m0_chirp_list*bin_bw + [-.5 .5]'*bin_bw;
  else
    f_chirp_offset = m0_chirp_list*bin_bw + f_offset_list';
  end
  f_chirp_range = max(max(f_chirp_offset)) - min(min(f_chirp_offset));
  f_total = n_chirp*f_incr_offset + f_chirp_range;
  f_range = max(freq) - min(freq);

  f1_start = -f_total/2;
  [temp,i_freq] = min(abs(freq-f1_start));
  f1_start = freq(i_freq);

  i_chirp = 0;
  m0_chirp_all = [];
  f_offset_all = [];
  df_dt_all = [];
  f1_all = [];

  for m0_chirp = m0_chirp_list
    for m0_offset = m0_offset_list
      for f_offset = f_offset_list

        i_chirp = i_chirp + 1;

        if isnan(m0_offset)
          m0_chirp_all(i_chirp) = m0_chirp + (rand(1)-.5);
        else
          m0_chirp_all(i_chirp) = m0_chirp + m0_offset;
        end
          
        if isnan(f_offset)
          f_offset_all(i_chirp) = (rand(1)-.5)*bin_bw;
        else
          f_offset_all(i_chirp) = f_offset;
        end

        df_dt_all(i_chirp) = m0_chirp_all(i_chirp)*bin_bw/t_max;
        f1_all(i_chirp) = f1_start + (i_chirp-1)*f_incr_offset + f_offset_all(i_chirp);
      end
    end
  end

  f2_all = f1_all + df_dt_all*t_max;
  
  fprintf('n_chirp=%.0f total range %.0f to %.0f vs %.0f to %.0f\n',...
        n_chirp,min(f1_all),max(f2_all),min(freq),max(freq));

  if (f_total > f_range)
    error(sprintf(...
      'Error in gen_multi_chirp_parameters, total freq range=%.0f > PFB range = %.0f ',...
      f_total,f_range));
  end

end
