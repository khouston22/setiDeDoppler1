function ddt_db = delta_dt_db(N_avg);
%
% function to add to the detection threshold for non-coherent summation
% of Stokes I powers
% for values of N_avg below 100
% returns 0 for N_avg above 100
%

ddt_db = zeros(size(N_avg));

ii = find(N_avg<100);

ddt_db(ii)  = interp1(log10([1 10 100]),[4.1 1.1 0],log10(N_avg(ii)),'spline');

