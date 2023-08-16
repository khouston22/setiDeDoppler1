function text_sc(xx,yy,text_str,font_size)
% function text_sc(xx,yy,text_str,font_size)
%
% writes text on current plot in screen coordinates
%
% x and y are between 0 and 1
%
% font_size is optional 
%

v = axis;
x = v(1) + xx*(v(2)-v(1));
y = v(3) + yy*(v(4)-v(3));

if (~exist('font_size'))
  text(x,y,text_str);		
else
  text(x,y,text_str,'FontSize',font_size);
end;
