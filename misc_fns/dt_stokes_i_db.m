function dt_db = dt_stokes_i_db(N_avg);
%
% function for the detection threshold for non-coherent summation N_avg times (=BT)
% (Stokes I = magnitude sum of 2 polarizations) 
% Approximates Chi Square for 4 N_avg degrees of freedom, Pd=.5 Pfa=1e-12
%

dt_db = 8 - 5*log10(N_avg) + delta_dt_db(N_avg);


