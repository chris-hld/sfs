function fontname(fname)
%FONTSIZE sets the figures font
%
%   Usage: fontname(fname)
%
%   Input options:
%       fname   - font name
%
%   FONTNAME(fname) sets the font of the active figure to the given type.
%
%   See also: GraphDefaults, fontsize

% FIXME: check together with GraphDefaults if the function is needed
% anymore

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************

% AUTHOR: Sascha Spors, Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));

if ~ischar(fname)
    error('%s: fname has to be a string!',upper(mfilename));
end


%% ===== Apply settings ==================================================
% Get handle for active figure
h=gca;
% Set the font type
set(h,'FontName',fname);
temp=get(h,'Xlabel');
xlabel(get(temp,'String'));
set(temp,'FontName',fname);
temp=get(h,'Ylabel');
ylabel(get(temp,'String'));
set(temp,'FontName',fname);
temp=get(h,'Title');
title(get(temp,'String'));
set(temp,'FontName',fname);
