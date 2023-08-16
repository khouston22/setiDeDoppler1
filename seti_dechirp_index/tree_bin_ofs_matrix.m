function [tree_bin_ofs_N,m_list] = tree_bin_ofs_matrix(N,N0,Lf,pos_neg)
%
% Function to compute turbo seti deDoppler bin offset matrix using fast
% Taylor tree algorithm
% 
% Outputs
%
% inputs
%
% N               1 x 1         last  stage number of time samples (power of 2)
% N0              1 x 1         first stage number of time samples (power of 2)
% Lf              1x1           =1 for standard critically sampled FB ("1x")
%                               =2 for 50% overlapped frequency bins ("2x")
%                               =4 for 50% overlapped frequency bins ("4x")
% pos_neg         1 x 1         =1 for positive drift rates (m=0..n_time-1)
%                               =-1 for negative drift rates (m=-n_time+1..0)
%                               =0 for pos & neg drift rates (m=-n_time+1..n_time-1) (default)
%
% outputs
%
% tree_bin_ofs_N  N x N          (m,n) or (drift index,time_index) 
%                                shift in bins at 1x as a function of t=nT
%                                required to compensate for drift rate m
%

  if (~exist('Lf','var')),  Lf = 1; end
  if isempty(Lf),           Lf = 1; end

  if (~exist('pos_neg','var')),  pos_neg = 0; end
  if isempty(pos_neg),           pos_neg = 0; end

  if N~=2.^round(log2(N))
    error(sprintf('Error in tree_bin_ofs_matrix, N not a power of 2, %.0f\n',N));
  end

  if N0~=2.^round(log2(N0))
    error(sprintf('Error in tree_bin_ofs_matrix, N0 not a power of 2, %.0f\n',N0));
  end

  for N1 = 2.^[log2(N0):log2(N)]
    
    freq_bin_offset_N1 = [0:N1-1]'*[0:N1-1]/(N1-1);
    int_bin_offset_N1  = round(freq_bin_offset_N1);
    
    if (N1==N0)
      if (Lf==1)
        tree_N1 = int_bin_offset_N1;
      else
        tree_N1 = round(Lf*freq_bin_offset_N1)/Lf;
      end
    else
      tree_N1_2 = tree_N1;
      tree_N1 = zeros(N1,N1);
      for i=1:N1/2
        tree_N1(2*i-1,1:N1/2) = tree_N1_2(i,:);
        tree_N1(2*i,1:N1/2)   = tree_N1_2(i,:);
        ii = N1/2+1:N1;
        tree_N1(2*i-1,ii) = tree_N1_2(i,:)+i-1;
        tree_N1(2*i,ii)   = tree_N1_2(i,:)+i;
      end
    end
  end
  
  if (pos_neg==1)
    m_list  = [0:N-1];
    tree_bin_ofs_N = tree_N1;
  elseif (pos_neg==-1)
    m_list  = [-(N-1):0];
    tree_bin_ofs_N = -tree_N1;
  else
    m_list  = [-(N-1):N-1];
    tree_bin_ofs_N = [-tree_N1(:,N:-1:2) tree_N1];
  end

end

