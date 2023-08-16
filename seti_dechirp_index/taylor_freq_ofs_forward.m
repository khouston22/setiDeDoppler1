function [tree_bin_ofs_N,rms_error_N,error_N] = ...
                             taylor_freq_ofs_forward(Nt,Nt0,Nr0,Lf,tree_bin_ofs0)
  %
  % calculates frequency offset matrix for each time step and freq rate
  % for Taylor De-Doppler algorithm
  % applies Taylor algorithm with initial stage N0 with unique
  % initial stage delays
  %
  % inputs either 
  % - a single Nr0 x Nt0 offset matrix (all input groups have identical offsets)
  % - an Nr0 x N offset matrix (each input group has unique offsets)
  %
  
  %
  % fill out tree_bin_ofs_N0 matrix to be N0 x N
  %
  
  n_stage = round(log2(Nt/Nt0));

  if Nt0*(2.^n_stage)~=Nt
    error(sprintf(...
      'Error in taylor_freq_ofs_forward, Nt/Nt0 not a power of 2, Nt=%.2f Nt0=%.0f\n',...
      Nt,Nt0));
  end
  
  [Nrow,Ncol] = size(tree_bin_ofs0);
  tree_bin_ofs0 = repmat(tree_bin_ofs0,1,Nt/Ncol);
  
  if (Nt==Nt0)
    tree_bin_ofs_N = tree_bin_ofs0;
    error_N = tree_bin_ofs_N - freq_bin_ofs(Nt0,Nr0);
    rms_error_N = rms(error_N(:));
    return
  end
  
  Nr = Nr0*(2.^n_stage);
 
  tree_bin_ofs_Ni2 = tree_bin_ofs0;
  
  for i_stage = 1:n_stage
    Nti = (2.^i_stage)*Nt0;
    Nri = (2.^i_stage)*Nr0;
    Nti2 = Nti/2;
    Nri2 = Nri/2;
    %dstageN = [[zeros(Nri2,Nti2);zeros(Nri2,Nti2)] [0:Nri2-1 1:Nri2]'*ones(1,Nti2)*Nt0/Nr0];
    %dstageN = [[zeros(Nri2,Nti2);zeros(Nri2,Nti2)] round([0:Nri2-1 1:Nri2]'*ones(1,Nti2)*Nt0/Nr0)];
    dstageN = [[zeros(Nri2,Nti2);zeros(Nri2,Nti2)] round([0:Nri2-1 1:Nri2]'*ones(1,Nti2)*Nt0/Nr0*Lf)/Lf];

    i12 = [1:2:Nri 2:2:Nri];
    [temp,i21] = sort(i12);
    P12 = eye(Nri); 
    P12 = P12(i12,:); % normal to [even rows ; odd rows]
    P21 = eye(Nri); 
    P21 = P21(i21,:); % [even rows ; odd rows] to normal
    
    tree_bin_ofs_Ni = zeros(Nri,Nti);
    
    n_grid = Nt/Nti;
    for j1 = 1:n_grid
      jj = (j1-1)*Nti + [1:Nti];
      QQ = tree_bin_ofs_Ni2(:,jj);
      % forward process
      tree_bin_ofs_Ni(:,jj) = P21*([eye(Nri2);eye(Nri2)]*QQ + dstageN);
    end
    tree_bin_ofs_Ni2 = tree_bin_ofs_Ni;
  end
  
  tree_bin_ofs_N = tree_bin_ofs_Ni;
  
  error_N = tree_bin_ofs_N - freq_bin_ofs(Nt,Nr);
  rms_error_N = rms(error_N(:));
 
  return
end
