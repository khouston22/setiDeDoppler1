function cbar_axes = colorbar_label(colorbar_location,label_string);
%function cbar_axes = colorbar_label(colorbar_location,label_string);
%
% creates colorbar on current plot, with text label
%
% inputs
%
% colorbar_location	string		location per colorbar function
%					e.g. 'East', 'EastOutside', etc.
% label_string		string		text for label
%
% outputs
%
% cbar_axes		handle		handle to axes to permit further manipulation
%

if (~exist('colorbar_location')),  colorbar_location='EastOutside'; end;
if (length(colorbar_location)==0), colorbar_location='EastOutside'; end;
if (~exist('label_string')),  label_string=''; end;
if (length(label_string)==0), label_string=''; end;

cbar_axes = colorbar(colorbar_location);
if (length(label_string)>0)

  if (strmatch('east',lower(colorbar_location)))

    set(get(cbar_axes,'YLabel'),'String',label_string);

  elseif (strmatch('west',lower(colorbar_location)))

    set(get(cbar_axes,'YLabel'),'String',label_string);

  elseif (strmatch('vertical',lower(colorbar_location)))

    set(get(cbar_axes,'YLabel'),'String',label_string);

  elseif (strmatch('north',lower(colorbar_location)))

    set(get(cbar_axes,'XLabel'),'String',label_string);

  elseif (strmatch('south',lower(colorbar_location)))

    set(get(cbar_axes,'XLabel'),'String',label_string);

  elseif (strmatch('horiz',lower(colorbar_location)))

    set(get(cbar_axes,'XLabel'),'String',label_string);

  end;
end;