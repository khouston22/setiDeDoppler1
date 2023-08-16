function out_string = get_path(in_string)
% function out_string = get_path(in_string)
%
% function to return directory from a full path\file name
% excludes final \ or /
%

out_string = [];

if (length(in_string) == 0), return; end;

ii = find((abs(in_string)==abs('\'))|(abs(in_string)==abs('/')));

if (length(ii)>0)
  out_string = in_string(1:max(ii)-1);
end;

