function [x_out,fs_out,f_bin,t_out,xcoef] = run_pfb_Lfx(x_in,coef,fs_in,n_out,pLf,qLf);
%
% function to run Polyphase Filter Bank with input x_in
% Generates interpolated output bins
%
% inputs
%
% x_in         n_in x 1   input waveform
% coef         MxN_tap    coefficients
% fs_in        1x1        PFB input sample rate, Hz
% n_out        1x1        number of PFB output samples per subfilter
% pLf,qLf      1x1        Lf = pLf/qLf filter bank overlap factor
%                         =1=1/1 for standard critically sampled FB ("1x")
%                         =2=2/1 for 50% overlapped frequency bins ("2x")
%                         rational number: Lf = pLf/qLf, e.g. Lf=3/2=1.5
% outputs
%
% x_out        n_freq x n_out  output waveform
% fs_out       1x1             filter bank output sample rate, Hz = fs_in/M
% f_bin        n_freq x 1      bin center freq values
% t_out        n_out x 1       output time adjusted for PFB delay N_tap/2/fs_out
% xcoef        M x N_tap x pLf extended coeff array (for debug)
%

Lf = pLf/qLf;

[M,N_tap] = size(coef);
n_freq = M*Lf;
N_fft = M/qLf;

if mod(n_freq,1)~=0
  error(sprintf('Error in run_pfb_Lfx, n_freq=M*Lf must be integer, %.1f=%.0f*%.2f\n',n_freq,M,Lf));
end

x_in = x_in(:);
n_in = length(x_in);
n_in_min = n_out*M + (N_tap-1)*M;

if (n_in < n_in_min)
  x_in = [x_in ; zeros(n_in_min-n_in,1)];
else
  x_in = x_in(1:n_in_min);
end

x_in = reshape(x_in,M,[]);
x_out = complex(NaN(n_freq,n_out,class(x_in)));

xcoef = complex(zeros(M,N_tap,pLf,class(x_in)));

for i1 = 0:pLf-1
  xcoef(:,:,i1+1) = (coef/M).*(exp(-1j*2*pi*i1*[0:M-1]'/n_freq)*exp(-1j*2*pi*[0:N_tap-1]*i1/Lf));
end

x_out1 = zeros(pLf,N_fft,class(x_in));
for i_out = 1:n_out
  xx = x_in(:,i_out:i_out+N_tap-1);
  for i1 = 0:pLf-1
    xxc = xx.*xcoef(:,:,i1+1);
    S1 = sum(xxc,2);
    S2 = sum(reshape(S1,N_fft,[]),2);
    x_out1(i1+1,:) = fftshift(fft(S2.',N_fft));
  end
  x_out(:,i_out) = x_out1(:);
end

fs_out = fs_in/M;
f_bin = [-Lf*M/2:Lf*M/2-1]'/(Lf*M)*fs_in;
t_out = ([0:n_out-1]+N_tap/2)/fs_out;

