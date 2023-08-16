function [xx_avg,t_avg] = pre_DD_avg(xx,t_pfb,N_pre_DD);
%
% Function to reduce spectrograms in the time dimension by averaging N_pre_DD
% successive lines in the spectrogram
% 
% Outputs
%
% inputs
%
% xx           n_freq1 x n_time  mag squared time waveforms for each frequency bin (power)
% t_pfb        n_time x 1        time sec for time series in each pfb bin
% N_pre_DD     1 x 1             number of spectrom time lines to average
%
% outputs
%
% xx_avg       n_freq1 x n_time2  deDoppler-shifted spectrogram (power)
% t_avg        n_time x 1        time sec for averaged lines in each pfb bin
%

[n_freq1,n_time] = size(xx);

T_line = (t_pfb(2) - t_pfb(1))*N_pre_DD;

n_time2 = ceil(n_time/N_pre_DD);

if (n_time2*N_pre_DD > n_time)
  n_extra = n_time2*N_pre_DD - n_time;
  xx = [xx zeros(n_freq1,n_extra)];
end

xx_avg = zeros(n_freq1,n_time2);
t_avg = t_pfb(1) + T_line/2 + T_line*[0:n_time2-1]';

for i = 1:n_time2
  ii = (i-1)*N_pre_DD + [1:N_pre_DD];
  xx_avg(:,i) = mean(xx(:,ii),2);
end 





