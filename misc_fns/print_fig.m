function p = print_fig(p);
%
% prints current figure to png or jpg file
%
% inputs
%
% p                    struct 
% p.output_directory   string   output directory including end slash
% p.file_name          string   base file name
% p.file_ext           string   '.png', 'jpg', '.emf'
% p.file_count         1x1      current count for file name suffix, will be incremented
% p.dpi                1x1      dots per inch
%
% outputs
%
% p                    struct   with updated file_count
%
% example
%
% p = print_fig();     % get structure
% p.output_directory = './my_directory/';
% p.file_name = 'out_name';
% for i=1:n_fig
%   (create figure)
%   p = print_fig(p);  % creates file in ./my_directory/out_name-[1..n_fig].png
% end

if (~exist('p','var')),  p=[]; end;
if (length(p)==0)
  p = struct(...
    'output_directory','', ...    % output directory including end slash
    'file_name','', ...           % base file name
    'file_ext','.png', ...        % extension
    'file_count',0, ...           % count of figures in succession
    'dpi',150, ...                % dots per inch
    'jpgnn',90);                  % jpeg quality (ignored for other types)
  return;
end

if ~exist(p.output_directory,'dir'), mkdir(p.output_directory); end

p.file_count = p.file_count + 1;

if (p.file_name(end)=='-')
  full_file_name = sprintf('%s%s%.0f%s',p.output_directory,p.file_name,p.file_count,p.file_ext);
else
  full_file_name = sprintf('%s%s-%.0f%s',p.output_directory,p.file_name,p.file_count,p.file_ext);
end

dpi_string = sprintf('-r%.0f',p.dpi);

if strcmp(p.file_ext,'.png')
  try
    print('-dpng',dpi_string,full_file_name);
  catch
    % repeat if error occurs, sometimes problems opening file
    print('-dpng',dpi_string,full_file_name);
  end
elseif strcmp(p.file_ext,'.jpg')
  jpeg_string = sprintf('-djpeg%.0f',p.jpgnn);
  try
    print(jpeg_string,full_file_name);
  catch
    % repeat if error occurs, sometimes problems opening file
    print(jpeg_string,full_file_name);
  end
elseif strcmp(p.file_ext,'.emf')
  try
    print('-dmeta',dpi_string,full_file_name);
  catch
    % repeat if error occurs, sometimes problems opening file
    print('-dmeta',dpi_string,full_file_name);
  end
end
