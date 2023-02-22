function [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% CBREWER - This function produces a colorbrewer table (rgb data) for a 
% given type, name and number of colors of the colorbrewer tables. 
% For more information on 'colorbrewer', please visit
% http://colorbrewer2.org/
% 
% The tables were generated from an MS-Excel file provided on the website
% http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_updates.html
%
% 
% [colormap]=cbrewer(ctype, cname, ncol, interp_method)
%
% INPUT:
%   - ctype: type of color table 'seq' (sequential), 'div' (diverging), 'qual' (qualitative)
%   - cname: name of colortable. It changes depending on ctype.
%   - ncol:  number of color in the table. It changes according to ctype and
%            cname
%   - interp_method: interpolation method (see interp1.m). Default is "cubic" )
% 
% A note on the number of colors: Based on the original data, there is
% only a certain number of colors available for each type and name of
% colortable. When 'ncol' is larger then the maximum number of colors
% originally given, an interpolation routine is called (interp1) to produce 
% the "extended" colormaps.
%
% Example:  To produce a colortable CT of ncol X 3 entries (RGB) of 
%           sequential type and named 'Blues' with 8 colors:
%                   CT=cbrewer('seq', 'Blues', 8);
%           To use this colortable as colormap, simply call:
%                   colormap(CT)
% 
%           To see the various colormaps available according to their types and
%           names, simply call: cbrewer()
%
%  This product includes color specifications and designs developed by
%  Cynthia Brewer (http://colorbrewer.org/).
%
% Author: Charles Robert
% email: tannoudji@hotmail.com
% Date: 06.12.2011
% ------------------------------
% 18.09.2015  Minor fixes, fixed a bug where the 'spectral' color table did not appear in the preview


% load colorbrewer data
load('colorbrewer.mat')
% initialise the colormap is there are any problems
colormap=[];
if (~exist('interp_method', 'var'))
    interp_method='cubic';
end

% If no arguments
if (~exist('ctype', 'var') | ~exist('cname', 'var') | ~exist('ncol', 'var'))
    disp(' ')
    disp('[colormap] = cbrewer(ctype, cname, ncol [, interp_method])')
    disp(' ')
    disp('INPUT:')
    disp('  - ctype: type of color table *seq* (sequential), *div* (divergent), *qual* (qualitative)')
    disp('  - cname: name of colortable. It changes depending on ctype.')
    disp('  - ncol:  number of color in the table. It changes according to ctype and cname')
    disp('  - interp_method:  interpolation method  (see interp1.m). Default is "cubic" )')
    
    disp(' ')
    disp('Sequential tables:')
    z={'Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd',...
             'Purples','RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'Spectral'};
    disp(z')     
         
    disp('Divergent tables:')
    z={'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'};
    disp(z')
    
    disp(' ')
    disp('Qualitative tables:')
    %getfield(colorbrewer, 'qual')
    z={'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'};
    disp(z')

    plot_brewer_cmap
    return
end

% Verify that the input is appropriate
ctype_names={'div', 'seq', 'qual'};
if (~ismember(ctype,ctype_names))
    disp('ctype must be either: *div*, *seq* or *qual*')
    colormap=[];
    return
end

if (~isfield(colorbrewer.(ctype),cname))
    disp(['The name of the colortable of type *' ctype '* must be one of the following:'])
    getfield(colorbrewer, ctype)
    colormap=[];
    return
end

if (ncol>length(colorbrewer.(ctype).(cname)))
%     disp(' ')
%     disp('----------------------------------------------------------------------')
%     disp(['The maximum number of colors for table *' cname '* is ' num2str(length(colorbrewer.(ctype).(cname)))])
%     disp(['The new colormap will be extrapolated from these ' num2str(length(colorbrewer.(ctype).(cname))) ' values'])
%     disp('----------------------------------------------------------------------')
%     disp(' ')
    cbrew_init=colorbrewer.(ctype).(cname){length(colorbrewer.(ctype).(cname))};
    colormap=interpolate_cbrewer(cbrew_init, interp_method, ncol);
    colormap=colormap./255;
    return
end

if (isempty(colorbrewer.(ctype).(cname){ncol}))
    
    while(isempty(colorbrewer.(ctype).(cname){ncol}))
        ncol=ncol+1;
    end        
    disp(' ')
    disp('----------------------------------------------------------------------')
    disp(['The minimum number of colors for table *' cname '* is ' num2str(ncol)])
    disp('This minimum value shall be defined as ncol instead')
    disp('----------------------------------------------------------------------')
    disp(' ')
end

colormap=(colorbrewer.(ctype).(cname){ncol})./255;

end