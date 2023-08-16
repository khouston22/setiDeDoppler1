function [GWmean_db] = GWmean(x_db,df_dt_list,df_dt_mean,df_dt_sigma);
%
% Function to produce a Gaussian-pdf-weighted average gain in dB
% where the pdf is a Gaussian over drift rate
% 
% Outputs
%
% inputs
%
% x_db         n_rate x n_case   dB gain
% df_dt_list   n_rate x 1        det_DD drift rate range Hz/sec
% df_dt_mean   1 x 1             drift rate mean (default=0)
% df_dt_sigma  1 x 1             drift rate std dev 
%                                default=max(abs(df_dt_list))/1.96
%                                i.e. df_dt_list covers 2 sigma, or 95% tail
%                                area if mean is zero
%
% outputs
%
% GWmean_db    1 x n_case        Gaussian weighted mean gain in dB
%

if ~exist('df_dt_mean'), df_dt_mean = 0; end
if isempty(df_dt_mean),  df_dt_mean = 0; end

if ~exist('df_dt_sigma'), df_dt_sigma = max(abs(df_dt_list))/1.96; end
if isempty(df_dt_sigma),  df_dt_sigma = max(abs(df_dt_list))/1.96; end

df_dt_list = df_dt_list(:);

n_rate = size(df_dt_list,1);

[m,n] = size(x_db);
if (m==n_rate)
  n_case = n;
elseif (n==n_rate)
  x_db = x_db.';
  n_case = m;
else
  error('Error in GWmean, x_db and df_dt_list don''t match');
end

w = exp(-.5*((df_dt_list-df_dt_mean)/df_dt_sigma).^2);

GWmean_db = NaN(1,n_case);

% compute weighted sum, excluding NaN data points

for i_case = 1:n_case
  ii = find(~isnan(x_db(:,i_case)));
  if (length(ii)>0)
    GWmean_db(i_case) = sum(x_db(ii,i_case).*w(ii))/sum(w(ii));
  end
end

