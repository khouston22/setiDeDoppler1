function [results_text,line_list,file_list,line_number] = ...
			grep(search_string,file_spec,n_context_lines)
% function [results_text,line_list,file_list,line_number] = ...
%			grep(search_string,file_spec,n_context_lines)
%
% looks for instances of the search string in the files indicated by
% the file spec string
% 
% inputs
%
% search_string		string		string to be searched in files
% file_spec		string		file spec e.g. '*.m'		default='*.m'
% 			OR dir struct	files from Dir function, e.g. D=dir('*.m')
%			OR cell array	strings in cell array with path\file name
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
%			D = dir('*.m'); grep('xyz',D)
%			D = dir('*.m'); D2 = struct2cell(D); D3=D2(1,:); grep('xyz',D3)
%			
% Note: use grep_subdir to search all subdirectories
%
% K Houston Oct 2004
% updated May 2006 to support grep_subdir
%

if (~exist('file_spec')),  file_spec='*.m'; end;
if (length(file_spec)==0), file_spec='*.m'; end;
if (~exist('n_context_lines')),  n_context_lines=2; end;
if (length(n_context_lines)==0), n_context_lines=2; end;

if (isstruct(file_spec))
  D2 = struct2cell(file_spec); 
  input_file_list=D2(1,:);
elseif (iscell(file_spec))
  input_file_list = file_spec;
else
  D = dir(file_spec);
  D2 = struct2cell(D); 
  file_spec_path = get_path(file_spec);
  if (length(file_spec_path)>0)
    file_spec_path = [file_spec_path '\'];
  end;

  input_file_list = [];
  for i=1:size(D2,2)
    input_file_list{i} = [file_spec_path D2{1,i}];
  end;
end;

results_text = [];
line_list = [];
file_list = [];
line_number = [];

for i_file = 1:length(input_file_list)
  input_file_name = input_file_list{i_file};

  fp = fopen(input_file_name,'r');

  line_history = [];
  n_line = 0;
  
  n_hits = 0;
  extra_line_count = 0;

  while (1)
    line1 = fgetl(fp);
    if (~isstr(line1))
      break;
    end;
    n_line = n_line + 1;
    if (extra_line_count>0)
      results_text = [results_text sprintf('%s\n',line1)];
      extra_line_count = extra_line_count - 1;
    elseif (length(strfind(line1,search_string))>0)
      extra_line_count = n_context_lines;
      results_text = [results_text ...
			sprintf('\n*** file %s, line %1.0f: \n\n',input_file_name,n_line)];
      for i=1:length(line_history)
        results_text = [results_text sprintf('%s\n',line_history{i})];
      end;
      line_history = [];
      results_text = [results_text sprintf('%s\n',line1)];
      n_hits = n_hits + 1;
      line_list{length(line_list)+1,1} = line1;
      file_list{length(file_list)+1,1} = input_file_name;
      line_number = [line_number ; n_line];
    elseif (n_context_lines>0)
      line_history{length(line_history)+1,1} = line1;
      n_lh = length(line_history);
      if (n_lh>n_context_lines)
        line_history = line_history(2:n_lh);
      end;
    end;
  end;
  %fprintf('Searched file %s, %1.0f lines, %1.0f hits\n',input_file_name,n_line,n_hits);

  fclose(fp);
end;