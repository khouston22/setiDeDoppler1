function [x_in,t_max,t_in] = gen_chirp_wf(fs_in,n_in,f1,df_dt,phi,class_name);
%
% function to generate chirp that is input to Polyphase Filter Bank
%
% inputs
%
% fs_in        1x1        sample rate Hz
% n_in         1x1        number of samples
% f1           1x1        start frequency at t=0
% df_dt        1x1        drift rate Hz/sec
% phi          1x1        additional phase, radians
% class_name   string     output type, ='double' (default), 'single'
%
% outputs
%
% x_in         n_in x 1   waveform
% t_max        1 x 1      time duration of chirp sec
% t_in         n_in x 1   input time starting at t=0
%

if (~exist('phi')), phi = 0; end
if isempty(phi),    phi = 0; end

if (~exist('class_name','var')),  class_name = 'double'; end
if isempty(class_name),           class_name = 'double'; end

t_max = n_in/fs_in;

% perform calculations in double precision

t_in = [0:n_in-1]'/fs_in;

if (1)
  x_in = exp(1j*(2*pi*f1*t_in + pi*df_dt*t_in.^2));
else
  f2 = f1 + t_max*df_dt;
  f_center = (f1 + f2)*.5;
  delta_f = f2 - f1;
  x_in = exp(1j*(2*pi*f_center*t_in + pi*delta_f*(t_in.*(t_in/t_max - 1))));
end

if (phi~=0)
  x_in = x_in*exp(1j*phi);
end

if strcmp(class_name,'single')
  x_in = single(x_in);
  t_in = single(t_in);
end


