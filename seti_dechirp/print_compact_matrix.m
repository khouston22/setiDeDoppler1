function print_compact_matrix(name_string,A,n_left,n_right,fid);
%
% Function to print out a matrix in a compact form
%
% 
% Outputs
%
% inputs
%
% name_string  string      user-selected text message before matrix values     
% A            n x m       matrix to be printed
% n_left       1 x 1       number of digits to left of decimat point
% n_right      1 x 1       number of digits to right of decimat point
% fid          1 x 1       file id (optional), =1 to stdout
%
% examples
% print_compact_matrix('in',in,2,0); prints matrix in
%

  if (~exist('fid','var')),  fid = 1; end  % standard output
  if isempty(fid),           fid = 1; end


  fprintf(fid,'%s\n',name_string);
  [m,n]=size(A);
  if (n_right>0)
    p_spec = sprintf('%%%.0f.%.0ff ',n_left+n_right+1,n_right);
  else
    p_spec = sprintf('%%%.0f.0f ',n_left);
  end
  for i=1:m
    fprintf(fid,p_spec,A(i,:));
    fprintf(fid,'\n');
  end
end