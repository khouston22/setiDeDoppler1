function [tree_bin_ofs_N,rms_error_N,epct_N,error_N] = ...
                        taylor_freq_ofs_forward_ext0(N,N0,ext_bin_ofs_N0)
  %
  % calculates frequency offset matrix for each time step and freq rate
  % for Taylor De-Doppler algorithm
  % applies Taylor algorithm with initial stage N0 with unique
  % initial stage delays
  %
  % inputs a special initial NxN first-stage delay matrix
  % (more hypothetical interest than practical)
  %
  % assumes Lf=1 for now, N x N (zero and positive drift rates)
  %
  
  %
  % fill out extended matrix as necessary
  %
  
  [Nrow,Ncol] = size(ext_bin_ofs_N0);
  ext_bin_ofs_N0 = repmat(ext_bin_ofs_N0,N/Nrow,N/Ncol);
  
  if (N<=N0)
    tree_bin_ofs_N = ext_bin_ofs_N0;
    error_N = tree_bin_ofs_N - freq_bin_ofs0(N);
    rms_error_N = rms(error_N(:));
    epct_N = find_epct(error_N);
    return
  end
 
  n_stage = log2(N/N0);

  ext_bin_ofs_Ni = ext_bin_ofs_N0;
  
  for i_stage = 1:n_stage
    Ni = (2.^i_stage)*N0;
    Ni2 = Ni/2;
    
    n_grid = N/Ni;
    for i1 = 1:n_grid
      ii = (i1-1)*Ni + [1:Ni];
      for j1 = 1:n_grid
        jj = (j1-1)*Ni + [1:Ni];
        ext_QQ = ext_bin_ofs_Ni(ii,jj);
        [tree_bin_ofs_Ni] = taylor_forward_1stage_ext(Ni,ext_QQ);
        ext_bin_ofs_Ni(ii,jj) = tree_bin_ofs_Ni;
      end
    end
  end
  
  tree_bin_ofs_N = tree_bin_ofs_Ni;
  
  error_N = tree_bin_ofs_N - freq_bin_ofs0(N);
  rms_error_N = rms(error_N(:));
  epct_N = find_epct(error_N);
 
  return;
end

function [tree_bin_ofs_N] = taylor_forward_1stage_ext(N,ext_bin_ofs_N2)
 % assumes Lf=1 for now, N x N (zero and positive drift rates)
  
  N2 = N/2;

  dstageN = [[zeros(N2);zeros(N2)] [0:N2-1 1:N2]'*ones(1,N2)];
  P12 = eye(N); P12 = P12([1:2:N 2:2:N],:); % normal to [even rows ; odd rows]
  P21 = P12';
  % forward process
  tree_bin_ofs_N = P21*(ext_bin_ofs_N2 + dstageN);
    
  return
end

