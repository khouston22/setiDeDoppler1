function [x_mod,x_bits,bits,offset] = gen_simple_mod(fs_in,n_in,mod_params,class_name);
%
% function to generate random modulation of unit amplitude
% intention is to model spreading of spectrum by modulation
%
% inputs
%
% fs_in        1x1        sample rate Hz
% n_in         1x1        number of samples
% mod_params   struct     modulation parameters, [] for none
%   BPSK: mod_params.mod_type='BPSK',mod_params.T_sym=symbol period,
%      produces random +1/-1 bits with random offset from start.
%     The carrier is suppressed as in standard comms systems.
%   QPSK: mod_params.mod_type='BPSK',mod_params.T_sym=symbol period,
%      produces random (bI + 1j*bQ)/sqrt(2) bits with random offset from start.
%      where bI and bQ are a series of +1/-1 bits.
%     The carrier is suppressed as in standard comms systems.
%   C-BPSK1: mod_params.mod_type='C-BPSK1',mod_params.T_sym=symbol period,
%     mod_params.Am=modulation amplitude 0<Am<1,
%     produces random symbols with of amplitude 1+Am (+1) and 1-Am (0),
%     or a constellation [1 + Am(1) (1), 1 + Am(-1)(0)].  Modulation is in
%     the inphase direction, and is a form of amplitude modulation.
%     The carrier is not suppressed to aid in signal detection.
%   C-BPSK2: mod_params.mod_type='C-BPSK2',mod_params.T_sym=symbol period,
%     mod_params.Am=modulation amplitude 0<Am<1,
%     produces random symbols with quadrature constellation
%     [1 + j Am (1), 1 -j Am (0)]
%     The modulation is in quadrature from the carrier.
%   C-QPSK: mod_params.mod_type='C-QPSK',mod_params.T_sym=symbol period,
%   mod_params.Am=moduluation amplitude 0<Am<1,
%   produces random bits with 2 bits/symbol with constellation
%     [1 + (1+j)Am (11), 1 + (-1+j)Am(01), 1 + (-1-j)Am (00), 1 + (1-j)Am (01)]
%     This provides twice the bit rate from the BPSK case.
%
%   Also, can specify mod_params.T_start, mod_params.T_on, mod_params.T_off 
%   to control duty cycle
% class_name   string     output type, ='double' (default), 'single'
%
% outputs
%
% x_mod        n_in x 1   modulation waveform (complex) - at sample rate
% x_bits       n_in x 1   modulation bits (real or complex) - at sample rate
% bits         n_sym x 1  modulation bits (real or complex) - one per symbol
% offset       1 x 1      offset in samples to first symbol transition
%

if (~exist('class_name','var')),  class_name = 'double'; end
if isempty(class_name),           class_name = 'double'; end

x_mod = complex(ones(n_in,1,class_name),zeros(n_in,1,class_name));

if isempty(mod_params)
  return;
elseif ~isfield(mod_params,'mod_type')
  return;
elseif strcmp(upper(mod_params.mod_type),'NONE')
  return;
end

Ts = 1/fs_in;
samples_per_sym = round(mod_params.T_sym/Ts);
n_sym = ceil(n_in/samples_per_sym) + 1;
T_sym = samples_per_sym*Ts;

% empty symbols before Tx
n_sym_start = 0;
if isfield(mod_params,'T_start')
  if ~isempty(mod_params.T_start)
    n_sym_start = round(mod_params.T_start/T_sym);
  end
end

% active symbols for each Tx
n_sym_active = n_sym - n_sym_start;
if isfield(mod_params,'T_on')
  if ~isempty(mod_params.T_on)
    n_sym_active = round(mod_params.T_on/T_sym);
  end
end

% inactive symbols after each Tx
n_sym_off = n_sym - n_sym_start - n_sym_active;
if isfield(mod_params,'T_off')
  if ~isempty(mod_params.T_off)
    n_sym_off = min(n_sym_off,round(mod_params.T_off/T_sym));
  end
end

offset = floor(rand(1,1)*samples_per_sym);

if contains(upper(mod_params.mod_type),'QPSK')
  bits = 2*randi(2,2,n_sym,class_name)-3;      % pseudo random +1/-1 sequence
  bits = disable_bits(bits,n_sym_start,n_sym_active,n_sym_off);
  x_bits = ones(samples_per_sym,1,class_name)*(bits(1,:)+1j*bits(2,:));
else
  bits = 2*randi(2,1,n_sym,class_name)-3;      % pseudo random +1/-1 sequence
  bits = disable_bits(bits,n_sym_start,n_sym_active,n_sym_off);
  x_bits = ones(samples_per_sym,1,class_name)*bits;
end
x_bits = x_bits(:);
x_bits = x_bits(offset+[1:n_in]);
ii = find(abs(x_bits)>0);

x_mod = complex(zeros(n_in,1,class_name));

if strcmp(upper(mod_params.mod_type),'BPSK')
  x_mod(ii) = x_bits(ii);
elseif strcmp(upper(mod_params.mod_type),'QPSK')
  x_mod(ii) = x_bits(ii)/sqrt(2);
elseif strcmp(upper(mod_params.mod_type),'C-BPSK1')
  x_mod(ii) = 1 + mod_params.Am*x_bits(ii);
elseif strcmp(upper(mod_params.mod_type),'C-BPSK2')
  x_mod(ii) = 1 + 1j*mod_params.Am*x_bits(ii);
elseif strcmp(upper(mod_params.mod_type),'C-QPSK')
  x_mod(ii) = 1 + mod_params.Am*x_bits(ii);
end



function bits = disable_bits(bits,n_sym_start,n_sym_active,n_sym_off)
%
% function to form an on/off pattern (arbitary duty cycle) by disabling
% sections of a bit stream
% bit vector has +1/-1 values
% bits set to 0 if inactive
% repeating pattern of n_sym_active bits followed by n_sym_off inactive bits,
% with n_sym_start inactive bits at start
%

n_sym = size(bits,2);
if (n_sym_start>0)
  bits(:,1:n_sym_start) = 0;
end

i_sym = 1 + n_sym_start;
n_sym_remain = n_sym - n_sym_start;
while n_sym_remain>0
  i_sym = i_sym + n_sym_active;
  n_zero = max(0,min(n_sym_off,n_sym - i_sym + 1));
  if (n_zero<=0)
    break;
  end
  bits(:,i_sym:i_sym+n_zero-1) = 0;
  i_sym = i_sym + n_zero;
end




