function [xx_agc,bin_power_in,bin_power_out] = ...
          rfi_agc1(xx,N_avg_f,N_avg_t,weight_exp,bin_limit,bin_replace,bin_noise_magsq);
%
% Function to suppress interference in spectrograms by inverse scaling
% bins in spectrogram which exhibit low Doppler rate and/or impulsive
% behavior
% Effectively this applies "AGC"=automatic gain control, except power in 
% high-energy bins are suppressed. Weighting is based on a selectable
% exponent
% 
% Outputs
%
% inputs
%
% xx            n_freq x n_time  mag squared time waveforms for each frequency bin (power)
% N_avg_f       1 x 1            average length in freq 
%                                =1 no avg, =n_freq single avg over all freq,
%                                else moving average window length (assumed odd)
% N_avg_t       1 x 1            average length in time
%                                =1 no avg, =n_time single avg over all freq,
%                                else moving average window length (assumed odd)
% weight_exp    1 x 1            weighting exponent (pos), =2 baseline
%                                =2 for 1/sigma^2 weight, =4 for 1/sigma^4
% bin_limit     1 x 1            bin value for limiting operation, def=10
%                                set to zero if no limiting desired
% bin_replace   1 x 1            bin replacement value if over limit threshold
% bin_noise_magsq   1 x 1        noise-only mag square value for noise only
%                                without RFI, set to NaN to estimate
% Note: bin limit threshold   = bin_limit*bin_noise_magsq
%       bin replacement value = bin_replace*bin_noise_magsq
%
% outputs
%
% xx_agc        n_freq x n_time  power-suppressed spectrogram
% bin_power_in  n_freq x 1       avg spectrogram power before weighting
% bin_power_out n_freq x 1       avg spectrogram power after  weighting
%

[n_freq,n_time] = size(xx);

if (~exist('N_avg_f','var')),  N_avg_f = 1; end
if isempty(N_avg_f),           N_avg_f = 1; end

if (~exist('N_avg_t','var')),  N_avg_t = n_time; end
if isempty(N_avg_t),           N_avg_t = n_time; end

if (~exist('weight_exp','var')),  weight_exp = 2; end
if isempty(weight_exp),           weight_exp = 2; end

if (~exist('bin_limit','var')),  bin_limit = 10; end
if isempty(bin_limit),           bin_limit = 10; end

if (~exist('bin_replace','var')),  bin_replace = 2; end
if isempty(bin_replace),           bin_replace = 2; end

if (~exist('bin_noise_magsq','var')),  bin_noise_magsq = NaN; end
if isempty(bin_noise_magsq),           bin_noise_magsq = NaN; end

weight_exp = abs(weight_exp);

xx_agc = zeros(n_freq,n_time);

bin_power_in  = mean(xx,2);

%
% apply limit and replacement option
%

if (bin_limit>0)
  if isnan(bin_noise_magsq)
    bin_noise_magsq = 1.44*median(bin_power_in);
  end
  bin_limit_threshold = bin_limit*bin_noise_magsq;
  bin_replacement_value = bin_replace*bin_noise_magsq;
  
  xx(find(xx>bin_limit_threshold)) = bin_replacement_value;
else
  bin_limit_threshold = NaN;
  bin_replacement_value = NaN;
end

fprintf('\nIn rfi_agc1, magsq=%.1f, limit=%.1f, replace=%.1f\n',...
        bin_noise_magsq,bin_limit_threshold,bin_replacement_value);

%
% average over frequency dimension as specified
%

if (N_avg_f == 1)
  xx_avg = xx;
elseif (N_avg_f == n_freq)
  xx_avg = ones(n_freq,1)*mean(xx,1);
else
  xx_avg = moving_avg_freq(xx,N_avg_f);
end

%
% average over time dimension as specified
%

if (N_avg_t == 1)
  %xx_avg = xx_avg;
elseif (N_avg_t == n_time)
  xx_avg = mean(xx_avg,2)*ones(1,n_time);
else
  xx_avg = moving_avg_time(xx_avg,N_avg_t);
end

%
% weight input by averaged bins
%

if (weight_exp==2)
  xx_agc = xx./xx_avg;
else
  xx_agc = xx.*(xx_avg.^(-weight_exp/2));
end

bin_power_out = mean(xx_agc,1);


function xx_avg = moving_avg_time(xx,N_avg_t)
%
% compute moving average of power over each freq bin
% using cumsum method
%

[n_freq,n_time] = size(xx);

N_avg = 2*floor(N_avg_t/2)+1; % assure odd

if (N_avg > n_time)
  N_avg = n_time-1;
end

N_avg2 = floor(N_avg/2);

power_cumsum = cumsum(xx,2);

avg_power_lower = power_cumsum(:,N_avg)/N_avg;
avg_power_upper = (power_cumsum(:,n_time)-power_cumsum(:,n_time-N_avg))/N_avg;

i1 = 2+N_avg2;
i2 = n_time - N_avg2;

xx_avg = zeros(n_freq,n_time);

for i = 1:i1-1
  xx_avg(:,i) = avg_power_lower;
end 

for i = i1:i2
  xx_avg(:,i) = (power_cumsum(:,i+N_avg2)-power_cumsum(:,i-N_avg2-1))/N_avg;
end 

for i = i2+1:n_time
  xx_avg(:,i) = avg_power_upper;
end 


function xx_avg = moving_avg_freq(xx,N_avg_f)
%
% compute moving average of power over each freq bin
% using cumsum method
%

[n_freq,n_time] = size(xx);

N_avg = 2*floor(N_avg_f/2)+1; % assure odd

if (N_avg > n_freq)
  N_avg = n_freq-1;
end

N_avg2 = floor(N_avg/2);

power_cumsum = cumsum(xx,1);

avg_power_lower = power_cumsum(N_avg,:)/N_avg;
avg_power_upper = (power_cumsum(n_freq,:)-power_cumsum(n_freq-N_avg,:))/N_avg;

i1 = 2+N_avg2;
i2 = n_freq - N_avg2;

xx_avg = zeros(n_freq,n_time);

for i = 1:i1-1
  xx_avg(i,:) = avg_power_lower;
end 

for i = i1:i2
  xx_avg(i,:) = (power_cumsum(i+N_avg2,:)-power_cumsum(i-N_avg2-1,:))/N_avg;
end 

for i = i2+1:n_freq
  xx_avg(i,:) = avg_power_upper;
end 




  



