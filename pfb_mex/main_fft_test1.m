%
% script to verify C-code mex version of taylorDD3.c gives identical results to 
% matlab function taylorDD3.m 
%

clear;
%clear classes;

addpath('..\seti_dechirp');
addpath('..\misc_fns');
addpath('..\pfb');

test_LR = 0;        % =1  test fft_LR, =0 test kiss_fft
use_random = 1;     % =1  use random input data, =0 create drift lines
print_matrices = 0; % =1 print mex and .m outputs and differences

if (test_LR)
  % delete fft_LR_mex.mexw64
  fprintf('Testing fft_LR_mex.c\n');
  mex fft_LR_mex.c fft_LR.c
else
  % delete fft_LR_mex.mexw64
  fprintf('Testing kiss_fft_mex.c\n');
  mex kiss_fft_mex.c kiss_fft.c
end

%return;

if (1)
  n_trial = 10;
  n_col = 10;
  N_list = [512 1024 8192];
else
  n_trial = 1;
  n_col = 2;
  N_list = [64];
  print_matrices = 1;
end

for N = N_list
  for i_trial = 1:n_trial
  
    x_in = complex(round(10*randn(N,n_col,'single')),round(10*randn(N,n_col,'single')));

    % matlab lib function  

    t0 = tic;
    x_out_m = fft(x_in,N);
    t_elapsed_m = toc(t0);

    % mex taylor_DD
    t0 = tic;
    if (test_LR)
      x_out_mex = fft_LR_mex(x_in,N);
    else
      x_out_mex = kiss_fft_mex(x_in,N);
    end
    t_elapsed_mex = toc(t0);

    difference = x_out_mex - x_out_m;
    max_diff = max(max(abs(difference)));

    %if (print_matrices) || (max_diff>0)
    if (print_matrices)
      print_compact_matrix('real x_in',real(x_in.'),3,0);
      print_compact_matrix('imag x_in',imag(x_in.'),3,0);
      print_compact_matrix('real x_out_mex',real(x_out_mex.'),3,0);
      print_compact_matrix('imag x_out_mex',imag(x_out_mex.'),3,0);
      print_compact_matrix('real x_out_m',real(x_out_m.'),3,0);
      print_compact_matrix('imag x_out_m',imag(x_out_m.'),3,0);
      print_compact_matrix('difference',difference',3,0);
    end

    fprintf('N=%4.0f n-col=%3.0f max diff=%2.0f, %.3f mex vs %.3f m msec\n\n',...
      N,n_col,max_diff,t_elapsed_mex*1e3,t_elapsed_m*1e3);
  end
end


