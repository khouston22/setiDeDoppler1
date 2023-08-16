function coef = calc_sinc_window_coefs(M,N_tap,window_name,apply_sinc,class_name);
%
% function to generate Polyphase Filter Bank coefficients
% based on sinc(n)*window(n) format
%
% inputs
%
% M            1x1        Number of subfilters at 1x = FFT size
% N_tap        1x1        Taps per subfilter, =1 for FFT filter bank
% window_name  string     window type
% apply_sinc   1x1        =1 applies sinc.*window, 0 = window only
% class_name   string     output type, ='double' (default), 'single'
%
%
% outputs
%
% coef         MxN_tap        coefficients
%

if (~exist('window_name')), window_name = 'Hann'; end
if isempty(window_name),    window_name = 'Hann'; end

if (~exist('apply_sinc')), apply_sinc = 1; end
if isempty(apply_sinc),    apply_sinc = 1; end

if (~exist('class_name','var')),  class_name = 'double'; end
if isempty(class_name),           class_name = 'double'; end

window_name = lower(window_name);

if strcmp(window_name,'hann')
  wind = hann(N_tap*M);
elseif contains(window_name,'lpf')
	wind = design_dec_lpf(N_tap,M,window_name);
  apply_sinc = 0;
elseif strcmp(window_name,'hamming')
	wind = hamming(N_tap*M);
elseif strcmp(window_name,'kaiser')
	wind = kaiser(N_tap*M,2.5);
elseif strcmp(window_name,'blackman')
	wind = blackman(N_tap*M);
elseif strcmp(window_name,'blackmanharris')
	wind = blackmanharris(N_tap*M);
elseif strcmp(window_name,'nuttallwin')
	wind = nuttallwin(N_tap*M);
elseif strcmp(window_name,'flattopwin')
	wind = flattopwin(N_tap*M);
elseif strcmp(window_name,'chebwin')
	wind = chebwin(N_tap*M);
elseif strcmp(window_name,'tukey50')
	wind = tukeywin(N_tap*M,.5);
elseif strcmp(window_name,'taylorwin')
	wind = taylorwin(N_tap*M);
elseif strcmp(window_name,'gausswin')
	wind = gausswin(N_tap*M);
elseif strcmp(window_name,'bartlett')
	wind = bartlett(N_tap*M);
elseif strcmp(window_name,'barthannwin')
	wind = barthannwin(N_tap*M);
elseif strcmp(window_name,'triang')
	wind = triang(N_tap*M);
elseif strcmp(window_name,'parzenwin')
	wind = parzenwin(N_tap*M);
elseif strcmp(window_name,'bohmanwin')
	wind = bohmanwin(N_tap*M);
elseif strcmp(window_name,'rectwin')
	wind = rectwin(N_tap*M);
else
  error(sprintf('Error in calc_sinc_window_coefs, window %s not recognized\n',...
    window_name));
end

if apply_sinc
  coef = sinc([-N_tap/2*M:N_tap/2*M-1]'/M).*wind;
else
  coef = wind;
end

if (N_tap==1)
  if ~strcmp(lower(window_name),'rectwin')
    coef = coef*2;
  end
end

coef = reshape(coef,M,N_tap);

if strcmp(class_name,'single')
  coef = single(coef);
end

function wind = design_dec_lpf(N_tap,M,window_name)

% assumes window_name = 'lpfxxx'
% xxx passband freq factor x 100, eg 100=>nominal fp=fs_in/M
% rrr = passband p-p ripple dB x 100, e.g. 100=> 1 dB p-p, 50=>0.5 dB p-p
% ss = stopband atten dB, e.g. 25=> 25 dB stopband attenuation

fc = str2num(window_name(4:6))/100/M;

% [b,n_coef,pp_ripple_db1,reject_db1] = ...
% 	remez_lpf2(f_pass,f_stop,fs,pp_ripple_db,reject_db,L,...
% 			plot_enable,iter_count_max,verbose);
    
% [b_fir_cell,n_coef] = ML_filter_design_sinc_hamming(fs0,M,1,N_tap);

h = fir1(N_tap*M-1,fc);
%h = fir1(N_tap*M-1,fc,chebwin(N_tap*M,25));
%h = firls( N_tap*M-1, [0 2*fc 2*fc 1], [1 1 0 0]).*kaiser(N_tap*M,5)' ;
%h = firls( N_tap*M-1, [0 2*fc 2*fc 1], [1 1 0 0]);
wind = M*h'/sum(h);


