function [y,i,j] = min2d(A);
%function [y,i,j] = min2d(A);
%
% function to find min of an array and the indices where it occurs
%
% inputs
%
% A	m x n	input array
%
% outputs
%
% y	1x1	min(min(A))
% i	1x1	row index of global min
% j	1x1	column index of global min
%
% ie y = A(i,j)
%

[y1, i1] = min(A);

[y, j] = min(y1);

i = i1(j);

% swap indices if A is a row vector

if (size(A,1)==1)
  temp = i;
  i = j;
  j = temp;
end



