function gfmt = gspec(number,n_sig,max_digits,pre_text,post_text)
%
% returns a format string with a variable format analogous to %g
% "How %g format should work"
%
% string is used in fprintf FORMAT field
%
% n_sig is the number of significant digits
%
% if abs(number)<10.^(max_digits) or abs(number)<10.^-(max_digits-n_sig+1), 
% %e will be output instead
%
% idea is to prevent %f from having too many digits for a very large number
% or a very small number
%
% assumes number is single value 
%

if (~exist('n_sig','var')),        n_sig=3; end;
if isempty(n_sig),                 n_sig=3; end;
if (~exist('max_digits','var')),   max_digits=7; end;
if isempty(max_digits),            max_digits=7; end;
if (~exist('pre_text','var')),     pre_text=''; end;
if isempty(pre_text),              pre_text=''; end;
if (~exist('post_text','var')),    post_text=''; end;
if isempty(post_text),             post_text=''; end;

number1 = max(abs(number(:)));
if number1==0
  log_num = 0;
else
  log_num = floor(log10(number1));
end

if log_num > (max_digits-1) 
  
  prec = n_sig-1;
  gfmt = sprintf('%s%%.%.0fe%s',pre_text,prec,post_text);
  
elseif log_num < -(max_digits-n_sig-1)
  
  prec = n_sig-1;
  gfmt = sprintf('%s%%.%.0fe%s',pre_text,prec,post_text);
  
else
  
  prec = max(0,n_sig-log_num-1);
  gfmt = sprintf('%s%%.%.0ff%s',pre_text,prec,post_text);
  
end
  
