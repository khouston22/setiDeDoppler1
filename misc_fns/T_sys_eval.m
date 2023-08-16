function [T_sys,T_gal,T_rcvr,T_atm] = T_sys_eval(freq,T_rcvr0,percentile);
%
% system temperature model
% T_rcvr0 - receiver and cosmic background temp
% at center of microwave window, ~2 GHz
% uses model of Robert Braun et al,  “Anticipated Performance of the 
% Square Kilometre Array - Phase 1 (SKA1), Version 1.0”, 
% SKA Organisation, Jodrell Bank, Lower UK, 2019
%
% output is size n_freq x n_T_rcvr0 
%

if (~exist('percentile')), percentile = 50; end;
if (isempty(percentile)), percentile = 50; end;

freq = freq(:);
T_rcvr0 = T_rcvr0(:)';

f_GHz = freq*1e-9;
T_rcvr = ones(size(freq))*T_rcvr0;

if (percentile == 50)
	T_gal = 25.2*(f_GHz/.408).^-2.75;
	T_atm = .04*max(0,f_GHz-3).^2;
	T_sys = T_rcvr + T_gal + T_atm - .0778;
elseif (percentile == 10)
	T_gal = 17.1*(f_GHz/.408).^-2.75;
	T_atm = .04*max(0,f_GHz-3).^2;
	T_sys = T_rcvr + T_gal + T_atm - .0567;
elseif (percentile == 90)
	T_gal = 54.8*(f_GHz/.408).^-2.75;
	T_atm = .04*max(0,f_GHz-3).^2;
	T_sys = T_rcvr + T_gal + T_atm - .1428;
else
  error(sprintf('Error in T_sys_eval, percentile=%.0f not recognized\n',percentile));
end

