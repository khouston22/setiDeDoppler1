function [y,i,j] = max2d(A);
%function [y,i,j] = max2d(A);
%
% function to find max of an array and the indices where it occurs
%
% inputs
%
% A	m x n	input array
%
% outputs
%
% y	1x1	max(max(A))
% i	1x1	row index of global max
% j	1x1	column index of global max
%
% ie y = A(i,j)
%


[y1, i1] = max(A);

[y, j] = max(y1);

i = i1(j);

% swap indices if A is a row vector

if (size(A,1)==1)
  temp = i;
  i = j;
  j = temp;
end




