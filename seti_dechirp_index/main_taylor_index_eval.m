%
% script to check Taylor (TurboSETI) and fastDD recursion tree error characteristics
% Taylor delays are exact, high confidence
%

clear;

addpath('..\misc_fns');

output_directory=sprintf('./plots-taylor-fastDD-error/');
if ~exist(output_directory,'dir'), mkdir(output_directory); end

% %default_figure_position = [775 43 758 591];	% set figure size for png import to pptx
% default_figure_position = [1001 43 758 591];	% set figure size for png import to pptx
% 
% figure(1);
% set(1,'Position',default_figure_position);
% commandwindow;
% 
% png = print_fig;
% png.output_directory = output_directory;

%for Lf = [1 2]
Lf = [1];
pos_neg = 1;
N0 = 8;

%for N=[16 32 64 128 256 512 1024]
%for N=[16 32 64 128 256]
for N=[1024]
        
  N2 = N/2;
  fprintf('\nN = %.0f, N0 = %.0f\n',N,N0);

  freq_bin_ofs_N = freq_bin_ofs0(N);
  freq_bin_ofs_N2 = freq_bin_ofs0(N/2);
  freq_bin_ofs_N0 = freq_bin_ofs0(N0);
  int_bin_ofs_N  = round(freq_bin_ofs_N*Lf)/Lf;
  int_bin_ofs_N2 = round(freq_bin_ofs_N2*Lf)/Lf;
  int_bin_ofs_N0 = round(freq_bin_ofs_N0*Lf)/Lf;

  ext_nom_ofsR_N0 = repmat(freq_bin_ofs_N0,N/N0,N/N0);
  ext_nom_ofsQ_N0 = repmat(int_bin_ofs_N0,N/N0,N/N0);

  extc_nom_ofsR_N0 = repmat(freq_bin_ofs_N0,1,N/N0);
  extc_nom_ofsQ_N0 = repmat(int_bin_ofs_N0,1,N/N0);

  % quantized error without fast algorithm
  
  rmse_direct = rms(int_bin_ofs_N(:)-freq_bin_ofs_N(:));
  epct_direct = find_epct(int_bin_ofs_N(:)-freq_bin_ofs_N(:));

  % unquantized Taylor error with fast algorithm
  
  [tree_ofs_taylorR,rmse_taylorR,epct_taylorR,error_taylorR] = taylor_freq_ofs_forward0(N,N0,freq_bin_ofs_N0);

  % nominal Taylor N from int_bin_ofs_N0
  
  [tree_ofs_taylorQ,rmse_taylorQ,epct_taylorQ] = taylor_freq_ofs_forward0(N,N0,int_bin_ofs_N0);

  P12 = eye(N); P12 = P12([1:2:N 2:2:N],:); % normal to [even rows ; odd rows]
  P21 = P12';
  
  error_N0_Q = int_bin_ofs_N0 - freq_bin_ofs_N0;
  error_taylorR_even_odd = P12*error_taylorR;

  if (0)
    pos_neg = 1;
    [total_bin_ofs_Nx] = tree_bin_ofs_matrix(N,N0,Lf,pos_neg);
    [total_bin_ofs_Ny] = tree_bin_ofs_matrix2(N,N0,int_bin_ofs_N0);
    diff_ofs = tree_ofs_taylorQ - total_bin_ofs_Ny;
    fprintf('Error Check forward taylor N0 matrix %.0f\n',max(abs(diff_ofs(:))));
  end

  % reverse process 1

  [ext_inv0R] = taylor_reverse(N,N0,freq_bin_ofs_N);
  ext_inv0Q = round(ext_inv0R*Lf)/Lf;
  diff_ext_inv0Q = ext_inv0Q - ext_nom_ofsQ_N0;

  % verify ext_bin_ofs_N0 inverts freq_bin_ofs_N exactly 
  % (note that this solution is unrealizable)

  [tree_ofs_inv0R,rmse_inv0R,epct_inv0R] = taylor_freq_ofs_forward_ext0(N,N0,ext_inv0R);
  [tree_ofs_inv0Q,rmse_inv0Q,epct_inv0Q] = taylor_freq_ofs_forward_ext0(N,N0,ext_inv0R);

  if (1)
    % block average rows in first stage and assess error
    [extc_inv1R] = block_average(ext_inv0R,N0,N);
  else
    [extc_inv1RR] = block_average(ext_inv0R,N,N0);
    [extc_inv1R] = block_transpose(extc_inv1RR,N0,N0);
  end
  [tree_ofs_inv1R,rmse_inv1R,epct_inv1R] = taylor_freq_ofs_forward0(N,N0,extc_inv1R);
  
  % quantize block average and assess error
  extc_inv1Q = round(extc_inv1R*Lf)/Lf;
  diff_extc_inv1Q = extc_inv1Q - extc_nom_ofsQ_N0;

  [tree_ofs_inv1Q,rmse_inv1Q,epct_inv1Q] = taylor_freq_ofs_forward0(N,N0,extc_inv1Q);
  
  % reverse process 2 - float
  
  [extc_inv2R] = taylor_inverse(N,N0,freq_bin_ofs_N);
  extc_inv2Q = round(extc_inv2R*Lf)/Lf;
  diff_extc_inv2R = extc_inv2R - extc_nom_ofsR_N0;
  diff_extc_inv2Q = extc_inv2Q - extc_nom_ofsQ_N0;

  % check how well extc_inv2R and extc_inv2Q inverts freq_bin_ofs_N 

  [tree_ofs_inv2R,rmse_inv2R,epct_inv2R,error_inv2R] = taylor_freq_ofs_forward0(N,N0,extc_inv2R);
  [tree_ofs_inv2Q,rmse_inv2Q,epct_inv2Q] = taylor_freq_ofs_forward0(N,N0,extc_inv2Q);


  % svd-based processes
  % set up initial guess
  
  dKin0 = repmat(freq_bin_ofs_N0,1,N/N0);
  [U0,S00,V0] = svd(dKin0); 
  if (1)
    dKin0 = extc_inv2R;
    %dKin0 = extc_inv2R +repmat(-.25*eye(N0),1,N/N0);
    %dKin0 = extc_inv2R + .25*randn(N0,N);
  end
  if (0)
    dKin0 = repmat(taylor_freq_ofs_forward0(N0,N0/2,freq_bin_ofs0(N0/2)),1,N/N0);
  end
  if (0)
    temp = taylor_freq_ofs_forward0(N,N0,freq_bin_ofs0(N0));
    dKin0 = taylor_inverse(N,N0,temp);
  end
  if (1)
    % choose U0 V0 from initial estimate
    [U0,S00,V0] = svd(dKin0); 
  else
    % keep the U0 V0 from freq_bin_ofs_N0
    S00 = U0'*dKin0*V0; 
  end
  extc_svd0R = dKin0;
  extc_svd0Q = round(extc_svd0R*Lf)/Lf;
  [tree_ofs_svd0R,rmse_svd0R,epct_svd0R,error_svd0R] = taylor_freq_ofs_forward0(N,N0,extc_svd0R);
  [tree_ofs_svd0Q,rmse_svd0Q,epct_svd0Q] = taylor_freq_ofs_forward0(N,N0,extc_svd0Q);

  % svd 1: non-linear least squares

  ds = .1;
  n_iter = 10;
  s0 = diag(S00(:,1:N0));
  S0 = [diag(s0) zeros(N0,N-N0)];
  extc_svd1R = dKin0;
  error_svd1R = error_svd0R;

  for i_iter = 1:n_iter
    if (1)
      % update SVD each iteration
      [J0,U0,S0,V0,s0] = find_jacobian1(N,N0,extc_svd1R,freq_bin_ofs_N,ds);
    elseif (1)
      % keep U0, V0 each iteration
      [J0] = find_jacobian2(N,N0,U0,S0,V0,freq_bin_ofs_N,ds);
    end
    ds0 = -pinv(J0)*error_svd1R(:);
    s1 = s0 + ds0;
    S1 = [diag(s1) zeros(N0,N-N0)];
    extc_svd1R = U0*S1*V0';  
    s0 = s1;
    S0 = S1;
    [tree_ofs_svd1R,rmse_svd1R,epct_svd1R,error_svd1R] = taylor_freq_ofs_forward0(N,N0,extc_svd1R);
    %fprintf('Iter %2.0f rmse=%.3f rms(ds0)=%.3f\n',i_iter,rmse_svd1R,rms(ds0));
    %fprintf('ds0= '); fprintf('%.3f ',ds0); fprintf('\n');
  end
  extc_svd1Q = round(extc_svd1R*Lf)/Lf;
  [tree_ofs_svd1Q,rmse_svd1Q,epct_svd1Q] = taylor_freq_ofs_forward0(N,N0,extc_svd1Q);

  % svd 2: Nelder Meade

  s0 = diag(S00(:,1:N0));
  S0 = [diag(s0) zeros(N0,N-N0)];

  fun = @(s_in)objectivefcn1(s_in,U0,V0,N,N0);
  s1 = fminsearch(fun,s0);

  S1 = [diag(s1) zeros(N0,N-N0)];
  extc_svd2R = U0*S1*V0';  
  extc_svd2Q = round(extc_svd2R*Lf)/Lf;
  [tree_ofs_svd2R,rmse_svd2R,epct_svd2R,error_svd2R] = taylor_freq_ofs_forward0(N,N0,extc_svd2R);
  [tree_ofs_svd2Q,rmse_svd2Q,epct_svd2Q] = taylor_freq_ofs_forward0(N,N0,extc_svd2Q);

  % svd 3: random perturbations

  s0 = diag(S00(:,1:N0));
  S0 = [diag(s0) zeros(N0,N-N0)];
  for i=1:100
    s1 = s0 + 1.0*randn(N0,1);
    S1 = [diag(s1) zeros(N0,N-N0)];
    dKin1 = U0*S1*V0';  
    [tree_ofs_svdxR,rmse_svdxR,epct_svdxR] = taylor_freq_ofs_forward0(N,N0,dKin1);
    s1_all(:,i) = s1; 
    rmse_all(i) = rmse_svdxR;
  end
  [min_rmse,i_min] = min(rmse_all);
  S1 = [diag(s1_all(:,i_min)) zeros(N0,N-N0)];
  extc_svd3R = U0*S1*V0';  
  extc_svd3Q = round(extc_svd3R*Lf)/Lf;
  [tree_ofs_svd3R,rmse_svd3R,epct_svd3R,error_svd3R] = taylor_freq_ofs_forward0(N,N0,extc_svd3R);
  [tree_ofs_svd3Q,rmse_svd3Q,epct_svd3Q] = taylor_freq_ofs_forward0(N,N0,extc_svd3Q);

  if (0)
    if (N<=64)
      %print_compact_matrix('int_bin_ofs_N',int_bin_ofs_N,2,0);
      %print_compact_matrix('int_bin_ofs_N2',int_bin_ofs_N2,2,0);
      print_compact_matrix('int_bin_ofs_N0',int_bin_ofs_N0,2,0);
      print_compact_matrix('diff_extc_inv2R',diff_extc_inv2R,2,1);
    end
  end
    
%   t_start = tic;
%   for i=1:100
%     %[tree_ofs_inv2Q,rmse_inv2Q] = taylor_freq_ofs_forward0(N,N0,extc_inv2Q);
%     [temp] = tree_bin_ofs_matrix2(N,N0,extc_inv2Q(1:N0,1:N0));
%   end
%   t_total=toc(t_start)/100

  fprintf('\nN = %.0f, N0 = %.0f\n',N,N0); 
  fprintf('rms error Direct         Q %.3f,       Q %2.0f%%\n',rmse_direct,epct_direct);
  fprintf('rms error Taylor R %.3f Q %.3f, R %2.0f%% Q %2.0f%%\n',...
    rmse_taylorR,rmse_taylorQ,epct_taylorR,epct_taylorQ);
  fprintf('rms error Inv0   R %.3f Q %.3f, R %2.0f%% Q %2.0f%%\n',...
    rmse_inv0R,rmse_inv0Q,epct_inv0R,epct_inv0Q);
  fprintf('rms error Inv1   R %.3f Q %.3f, R %2.0f%% Q %2.0f%%\n',...
    rmse_inv1R,rmse_inv1Q,epct_inv1R,epct_inv1Q);
  fprintf('rms error Inv2   R %.3f Q %.3f, R %2.0f%% Q %2.0f%%\n',...
    rmse_inv2R,rmse_inv2Q,epct_inv2R,epct_inv2Q);
  fprintf('rms error svd0   R %.3f Q %.3f, R %2.0f%% Q %2.0f%%\n',...
    rmse_svd0R,rmse_svd0Q,epct_svd0R,epct_svd0Q);
  fprintf('rms error svd1   R %.3f Q %.3f, R %2.0f%% Q %2.0f%%\n',...
    rmse_svd1R,rmse_svd1Q,epct_svd1R,epct_svd1Q);
  fprintf('rms error svd2   R %.3f Q %.3f, R %2.0f%% Q %2.0f%%\n',...
    rmse_svd2R,rmse_svd2Q,epct_svd2R,epct_svd2Q);
  fprintf('rms error svd3   R %.3f Q %.3f, R %2.0f%% Q %2.0f%%\n',...
    rmse_svd3R,rmse_svd3Q,epct_svd3R,epct_svd3Q);

end

%         if (1)
%           png.file_name = sprintf('01-%s-error-hist-%s-N-%.0f-N0-%.0f',alg_str,Lf_str,N,N0);
%           png.file_count = 0;
%           figure(1); clf;
%           hist(bin_error_N(:),101);
%           v=axis; axis([-1 1 v(3:4)]);
%           xlabel('Fractional Error (Frequency Bins)');
%           ylabel('Relative Frequency');
%           title(sprintf('%s DeDoppler %s Fractional Error, N=%.0f, N0=%.0f',alg_str,Lf_str,N,N0));
%           text_sc(.05,.95,results_string);
%           grid;
%           png = print_fig(png);
%           %pause(2);
%         end
% 

function [bin_avg_Nx] = block_average(ext_bin_ofs_Nx,Nrow1,Ncol1)

  [Nrow,Ncol] = size(ext_bin_ofs_Nx);
  N_avg_r = Nrow/Nrow1;
  N_avg_c = Ncol/Ncol1;
  bin_avg_Nx = zeros(Nrow1,Ncol1);
  for i_avg_r = 1:N_avg_r
    ii = (i_avg_r-1)*Nrow1 + [1:Nrow1];
    for i_avg_c = 1:N_avg_c
      jj = (i_avg_c-1)*Ncol1 + [1:Ncol1];
      bin_avg_Nx = bin_avg_Nx + ext_bin_ofs_Nx(ii,jj); 
    end
  end
  bin_avg_Nx = bin_avg_Nx/N_avg_r/N_avg_c;
end

function [Abt] = block_transpose(A,Nrow1,Ncol1)

  [Nrow,Ncol] = size(A);
  N_blk_r = Nrow/Nrow1;
  N_blk_c = Ncol/Ncol1;
  Abt =zeros(Ncol,Nrow);
  for i_blk_r = 1:N_blk_r
    ii = (i_blk_r-1)*Nrow1 + [1:Nrow1];
    for i_blk_c = 1:N_blk_c
      jj = (i_blk_c-1)*Ncol1 + [1:Ncol1];
      Abt(jj,ii) = A(ii,jj); 
    end
  end
end

function epct = find_epct(error_N)

error_N = error_N(:);
ii = find(abs(error_N)>.5);

epct = length(ii)/length(error_N)*100;

end

function [tree_bin_ofs_N] = tree_bin_ofs_matrix2(N,N0,tree_bin_ofs_N0)
  % assumes Lf=1 for now, N x N (zero and positive drift rates)
  % applies Taylor algorithm with initial stage N0 with single common
  % initial stage delays
  % tree_bin_ofs_N0 should be N0 x N0 
  
  %freq_bin_ofs_N0 = repmat(freq_bin_ofs_N0,1,N/N0);

  n_stage = log2(N/N0);

  if (n_stage==0)
    tree_bin_ofs_N = tree_bin_ofs_N0;
    return;
  end

  tree_bin_ofs_N2 = tree_bin_ofs_N0;

  for i_stage = 1:n_stage
    N = (2.^i_stage)*N0;
    N2 = N/2;

    dstageN = [[zeros(N2);zeros(N2)] [0:N2-1 1:N2]'*ones(1,N2)];
    P12 = eye(N); P12 = P12([1:2:N 2:2:N],:); % normal to [even rows ; odd rows]
    P21 = P12';
    % forward process
    tree_bin_ofs_N = P21*([eye(N2);eye(N2)]*tree_bin_ofs_N2*[eye(N2) eye(N2)] + dstageN);
    
    tree_bin_ofs_N2 = tree_bin_ofs_N;
  end
  return;
end

function [J0,U0,S0,V0,s0] = find_jacobian1(N,N0,dKin0,tree_ofs_star,ds);

  [U0,S0,V0] = svd(dKin0); 
  s0 = diag(S0(:,1:N0));
  
  J0 = zeros(N^2,N0);
  
  for i = 1:N0
    s1 = s0;
    s1(i) = s1(i) + ds;
    dKin1 = U0*[diag(s1) zeros(N0,N-N0)]*V0';
    tree_ofs_svd = taylor_freq_ofs_forward0(N,N0,dKin1);
    J0(:,i) = (tree_ofs_svd(:)-tree_ofs_star(:))/ds;
  end

end

function [J0] = find_jacobian2(N,N0,U0,S0,V0,tree_ofs_star,ds);

  s0 = diag(S0(:,1:N0));
  
  J0 = zeros(N^2,N0);
  
  for i = 1:N0
    s1 = s0;
    s1(i) = s1(i) + ds;
    dKin1 = U0*[diag(s1) zeros(N0,N-N0)]*V0';
    tree_ofs_svd = taylor_freq_ofs_forward0(N,N0,dKin1);
    J0(:,i) = (tree_ofs_svd(:)-tree_ofs_star(:))/ds;
  end

end

function [ext_bin_ofs_N0] = taylor_reverse(N,N0,bin_ofs_N)
 % assumes Lf=1 for now, N x N (zero and positive drift rates)
  
  n_stage = log2(N/N0);

  if (n_stage==0)
    ext_bin_ofs_N0 = bin_ofs_N;
    return
  end

  [ext_bin_ofs_N0] = taylor_reverse_1stage(N,bin_ofs_N);
  
  for i_stage = n_stage-1:-1:1
    Ni = (2.^i_stage)*N0;
    Ni2 = Ni/2;
    
    n_grid = N/Ni;
    for i1 = 1:n_grid
      ii = (i1-1)*Ni + [1:Ni];
      for j1 = 1:n_grid
        jj = (j1-1)*Ni + [1:Ni];
        QQ = ext_bin_ofs_N0(ii,jj);
        [ext_QQ] = taylor_reverse_1stage(Ni,QQ);
        ext_bin_ofs_N0(ii,jj) = ext_QQ;
      end
    end
  end
  return
end

function [ext_bin_ofs_N2] = taylor_reverse_1stage(N,bin_ofs_N)
 % assumes Lf=1 for now, N x N (zero and positive drift rates)
  
  N2 = N/2;

  dstageN = [[zeros(N2);zeros(N2)] [0:N2-1 1:N2]'*ones(1,N2)];
  P12 = eye(N); P12 = P12([1:2:N 2:2:N],:); % normal to [even rows ; odd rows]
  P21 = P12';
  % reverse process
  ext_bin_ofs_N2 = (P12*bin_ofs_N - dstageN);
  
  return;
end



function [extc_bin_ofs_N0] = taylor_inverse(N,N0,bin_ofs_N)
 % assumes Lf=1 for now, N x N (zero and positive drift rates)
 % output is N0 x N
  
  n_stage = log2(N/N0);

  if (n_stage==0)
    extc_bin_ofs_N0 = bin_ofs_N;
    return
  end

  [extc_bin_ofs_Ni] = taylor_inverse_1stage(N,bin_ofs_N);
   
  if (n_stage==1)
    extc_bin_ofs_N0 = extc_bin_ofs_Ni;
    return
  end

  for i_stage = n_stage-1:-1:1
    Ni = (2.^i_stage)*N0;
    Ni2 = Ni/2;
    extc_bin_ofs_Nim1 = zeros(Ni2,N);
    n_grid = N/Ni/2;
    for j1 = 1:n_grid
      jj1 = (j1-1)*2*Ni + [1:Ni];
      jj2 = jj1 + Ni;
      Q1 = extc_bin_ofs_Ni(:,jj1);
      Q2 = extc_bin_ofs_Ni(:,jj2);
      [extc_Q1] = taylor_inverse_1stage(Ni,Q1);
      [extc_Q2] = taylor_inverse_1stage(Ni,Q2);
      extc_bin_ofs_Nim1(:,jj1) = extc_Q1;
      extc_bin_ofs_Nim1(:,jj2) = extc_Q2;
    end
    extc_bin_ofs_Ni = extc_bin_ofs_Nim1;
  end
  extc_bin_ofs_N0 = extc_bin_ofs_Ni;
  return
end


function [extc_bin_ofs_N2] = taylor_inverse_1stage(N,bin_ofs_N)
 % assumes Lf=1 for now, N x N (zero and positive drift rates)
 % output is N/2 x N
  
  N2 = N/2;

  dstageN = [[zeros(N2);zeros(N2)] [0:N2-1 1:N2]'*ones(1,N2)];
  P12 = eye(N); P12 = P12([1:2:N 2:2:N],:); % normal to [even rows ; odd rows]
  P21 = P12';
  
  % reverse process
  if (1)
    extc_bin_ofs_N2 = [eye(N2) eye(N2)]*(P12*bin_ofs_N - dstageN)/2;
  elseif (0)
    temp = (P12*bin_ofs_N - dstageN)*[eye(N2) ; eye(N2)]/2;
    extc_bin_ofs_N2 = block_transpose(temp,N2,N2);
  else
    taper2 = diag([0:N2-1])/N2;
    taper1 = eye(N2) - taper2;
    temp = (P12*bin_ofs_N - dstageN)*[taper1 ;taper2];
    extc_bin_ofs_N2 = block_transpose(temp,N2,N2);
  end
  return;
end

