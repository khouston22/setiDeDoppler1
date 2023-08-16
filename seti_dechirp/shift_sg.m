function [xx_ssg,f_ssg] = shift_sg(xx,freq1,t_out,f1,df_dt,df_extract);
%
% Function to shift a spectrogram xx in frequency over time to compensate for
% a drift rate df_dt - will shift to nearest neighbor
% 
% Outputs
%
% inputs
%
% f1           1 x 1             "truth" start frequency Hz
% df_dt        1 x 1             "truth" drift rate Hz/sec
% xx           n_freq1 x n_time  max squared time waveforms for each frequency bin (power)
% freq1        n_freq1 x 1       PFB bin center freq values
% t_out        n_time x 1        time sec for time series in each pfb bin
% df_extract   1 x 1             number of frequency bins for zoom (odd #, def=11)
%
% outputs
%
% xx_ssg       nfb1 x n_time  deDoppler-shifted spectrogram (power)
% f_ssg        nfb1 x 1          det_chirp bin center freq values
%
% where
%
% nfb1 = number of frequency bins for zoom (odd #, def=11)
%        2*round(df_extract/(df_bin1/2)+1
%

[n_freq1,n_time] = size(xx);
N = 2.^ceil(log2(n_time));

df_bin1 = freq1(2)-freq1(1);  % may be less than 1/Ts if bins are overlapped
Ts = t_out(2)-t_out(1);
fs = 1/Ts;
Lf = round(1/(df_bin1*Ts));

if (~exist('df_extract','var')),  df_extract = 10*fs; end
if isempty(df_extract),           df_extract = 10*fs; end

nfb2 = round(df_extract/df_bin1/2);
nfb1 = 2*nfb2 + 1;

df_extract = nfb1*df_bin1;

%
% 
%

f_chirp = f1 + df_dt*t_out;

if (max(f_chirp)+df_extract/2 >= max(freq1)) || ...
   (min(f_chirp)-df_extract/2 <= min(freq1))
  fprintf(...
'\nWarning in shift_sg:\nChirp %.1f-%.1f outside of spectrogram range, %.1f-%.1f\n\n',...
  min(f_chirp),max(f_chirp),min(freq1),max(freq1));
end

xx_ssg = zeros(nfb1,N);

for n = 1:n_time
  %[temp,i_freq] = min(abs(freq1-f_chirp(n)));
  i_freq = 1 + round((f_chirp(n)-freq1(1))/df_bin1);
  xx_ssg(:,n) = xx(i_freq+[-nfb2:nfb2],n);
end 

[temp,i_freq] = min(abs(freq1-f_chirp(1)));
f_ssg = freq1(i_freq+[-nfb2:nfb2]);





