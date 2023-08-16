function [results_text,line_list,file_list,line_number] = ...
			grep_subdir(search_string,file_spec,n_context_lines)
% function [results_text,line_list,file_list,line_number] = ...
%			grep_subdir(search_string,file_spec,n_context_lines)
%
% looks for instances of the search string in the files indicated by
% the file spec string
%
% looks in all subdirectories of current directory
% 
% inputs
%
% search_string		string		string to be searched in files
% file_spec		string		file spec e.g. '*.m'		default='*.m'
% 					note: dir struct and cell array not supported
% context_lines		1x1		number of lines before/after string match default=2
%
% outputs
%
% results_text		string		text output with file name and context lines
% line_list		N x 1 cell 	strings where search string found
% file_list		N x 1 cell 	file_names corresponding to line_list
% line_number		N x 1 		line number list corresponding to line_list
%
% simplest call: 	grep('xyz')	looks for string 'xyz' among *.m in current directory
% equivalents: 		grep('xyz','*.m')
%			
%
% K Houston May 2006
%

if (~exist('file_spec')),  file_spec='*.m'; end;
if (length(file_spec)==0), file_spec='*.m'; end;
if (~exist('n_context_lines')),  n_context_lines=2; end;
if (length(n_context_lines)==0), n_context_lines=2; end;

if (isstruct(file_spec))
  error('Error: dir struct file spec not supported for grep_subdir');
elseif (iscell(file_spec))
  error('Error: cell array file name list not supported for grep_subdir');
end;

file_spec_path = get_path(file_spec);
file_spec_file = get_file_name(file_spec);

if (length(file_spec_path)>0)
  D = dir(file_spec_path);	% directory from file spec
  file_spec_path = [file_spec_path '\'];
else
  D = dir;			% current directory
end;
D = D(3:end);			% eliminate '.','..' entries
D2 = struct2cell(D);

input_name_list = [];
for i=1:size(D2,2)
  input_name_list{i} = [file_spec_path D2{1,i}];
end;

input_isdir_list=D2(4,:);

%
% call grep for current directory
%

[results_text,line_list,file_list,line_number] = ...
			grep(search_string,file_spec,n_context_lines);

%
% call grep_subdir recursively for each subdirectory of current directory
%

for i_dir = 1:length(input_isdir_list)
  if (input_isdir_list{i_dir}==1)
    dir_name = input_name_list{i_dir};
    ext_file_spec = [dir_name '\' file_spec_file];

    %fprintf('file_spec_path=%s, ext_file_spec=%s\n',file_spec_path,ext_file_spec);
    fprintf('Checking %s\n',ext_file_spec);

    [results_text1,line_list1,file_list1,line_number1] = ...
		grep_subdir(search_string,ext_file_spec,n_context_lines);

    if (size(results_text1,1) > 0)
      results_text = char(results_text,results_text1);
    end;
    line_list    = [line_list ; line_list1];
    file_list    = [file_list ; file_list1];
    line_number  = [line_number ; line_number1];
  end;
end;

