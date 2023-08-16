function [tree_bin_ofs_N,rms_error_N,epct_N,error_N] = ...
                             taylor_freq_ofs_forward0(N,N0,tree_bin_ofs_N0)
  %
  % calculates frequency offset matrix for each time step and freq rate
  % for Taylor De-Doppler algorithm
  % applies Taylor algorithm with initial stage N0 with unique
  % initial stage delays
  %
  % inputs either 
  % - a single N0 x N0 offset matrix (all input groups have identical offsets)
  % - an N0 x N offset matrix (each input group has unique offsets)
  %
  % assumes Lf=1 for now, N x N (zero and positive drift rates)
  %
  
  %
  % fill out tree_bin_ofs_N0 matrix to be N0 x N
  %
  
  [Nrow,Ncol] = size(tree_bin_ofs_N0);
  tree_bin_ofs_N0 = repmat(tree_bin_ofs_N0,1,N/Ncol);
  
  if (N<=N0)
    tree_bin_ofs_N = tree_bin_ofs_N0;
    error_N = tree_bin_ofs_N - freq_bin_ofs0(N);
    rms_error_N = rms(error_N(:));
    epct_N = find_epct(error_N);
    return
  end
 
  n_stage = log2(N/N0);

  tree_bin_ofs_Ni2 = tree_bin_ofs_N0;
  
  for i_stage = 1:n_stage
    Ni = (2.^i_stage)*N0;
    Ni2 = Ni/2;
    dstageN = [[zeros(Ni2);zeros(Ni2)] [0:Ni2-1 1:Ni2]'*ones(1,Ni2)];
    P12 = eye(Ni); 
    P12 = P12([1:2:Ni 2:2:Ni],:); % normal to [even rows ; odd rows]
    P21 = P12';
    tree_bin_ofs_Ni = zeros(Ni);
    
    n_grid = N/Ni;
    for j1 = 1:n_grid
      jj = (j1-1)*Ni + [1:Ni];
      QQ = tree_bin_ofs_Ni2(:,jj);
      % forward process
      tree_bin_ofs_Ni(:,jj) = P21*([eye(Ni2);eye(Ni2)]*QQ + dstageN);
    end
    tree_bin_ofs_Ni2 = tree_bin_ofs_Ni;
  end
  
  tree_bin_ofs_N = tree_bin_ofs_Ni;
  
  error_N = tree_bin_ofs_N - freq_bin_ofs0(N);
  rms_error_N = rms(error_N(:));
  epct_N = find_epct(error_N);
 
  return
end
