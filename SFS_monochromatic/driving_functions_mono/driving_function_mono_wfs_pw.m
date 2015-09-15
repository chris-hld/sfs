function D = driving_function_mono_wfs_pw(x0,nx0,nk,f,conf)
%DRIVING_FUNCTION_MONO_WFS_PW returns the driving signal D for a plane wave in
%WFS
%
%   Usage: D = driving_function_mono_wfs_pw(x0,nx0,nk,f,[conf])
%
%   Input parameters:
%       x0          - position of the secondary sources / m [nx3]
%       nx0         - directions of the secondary sources / m [nx3]
%       nk          - direction of plane wave / m [nx3]
%       f           - frequency of the monochromatic source / Hz
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_WFS_PW(x0,nx0,nk,f,conf) returns WFS driving signals
%   for the given secondary sources, the virtual plane wave direction and the
%   frequency f.
%
%   References:
%       H. Wierstorf, J. Ahrens, F. Winter, F. Schultz, S. Spors (2015) -
%       "Theory of Sound Field Synthesis"
%
%   See also: driving_function_mono_wfs, driving_function_imp_wfs_ps

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
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
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 5;
narginchk(nargmin,nargmax);
isargmatrix(x0,nx0,nk);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
xref = conf.xref;
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% Frequency
omega = 2*pi*f;


if strcmp('2D',dimension) || strcmp('3D',dimension)

    % === 2- or 3-Dimensional ============================================

    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % D using a plane wave as source model
        %
        %              i w
        % D(x0,w) =  2 --- nk nx0  e^(-i w/c nk x0)
        %               c
        %
        % see Wierstorf et al. (2015), eq.(#D:wfs:pw)
        %
        D = -2*1i*omega/c .* vector_product(nk,nx0,2) .*  ...
            exp(-1i*omega/c.*vector_product(nk,x0,2));
        %
    elseif strcmp('delft1988',driving_functions)
        % --- Delft 1988 -------------------------------------------------
        to_be_implemented;
        %
    else
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a plane wave.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

  % === 2.5-Dimensional ================================================
  
  % Reference point
  xref = repmat(xref,[size(x0,1) 1]);
  switch driving_functions
    case {'default', 'reference_point'}
      % Driving function with only one stationary phase approximation, i.e.
      % reference to one point in field
      %        ______________
      % g0 = \| 2pi |xref-x0|
      %
      g0 = sqrt( 2*pi*vector_norm(xref-x0,2) );
      %
      % D_2.5D using a plane wave as source model
      %                             ___
      %                            |i w|
      % D_2.5D(x0,w) = g0 nk nx0 _ |---  e^(-i w/c nk x0)
      %                           \| c
      %
      D = 2*g0 .* vector_product(nk,nx0,2) .* sqrt(1i*omega/c) .* ...
        exp(-1i*omega/c.*vector_product(nk,x0,2));
    case {'reference_line'}
      % Driving function with two stationary phase approximations,
      % reference to a line parallel to a LINEAR secondary source distribution
      %
      % distance ref-line to linear ssd
      dref = vector_product(xref-x0,nx0,2);
      % scalar product of plane wave and ssd normal ( npw * nx0 )
      cospw = vector_product(nk,nx0,2);
      %
      % 2.5D correction factor
      %        _____________________________
      % g0 = \| 2pi * d_ref * ( npw * nx0 )
      g0 = sqrt( 2*pi*dref.*cospw );
      %
      % D_2.5D using a plane wave as source model
      %                        ___
      %                       |i w|
      % D_2.5D(x0,w) = 2* g0 _|---  e^(-i w/c nk x0)
      %                      \| c
      %
      D = 2*g0 .* sqrt(1i*omega/c) .* exp( -1i*omega/c.*cospw );
    case {'legacy'}
      % --- SFS Toolbox ------------------------------------------------
      % 2.5D correction factor
      %        ______________
      % g0 = \| 2pi |xref-x0|
      %
      g0 = sqrt(2*pi*vector_norm(xref-x0,2));
      %
      % D_2.5D using a plane wave as source model
      %                               ___
      %                              | w |
      % D_2.5D(x0,w) = 2 g0 nk nx0 _ |---  e^(-i w/c nk x0)
      %                             \|i c
      %
      % see Wierstorf et al. (2015), eq.(#D:wfs:pw:2.5D)
      %
      D = -2*g0 .* vector_product(nk,nx0,2) .* sqrt(omega/(1i*c)) .* ...
        exp(-1i*omega/c.*vector_product(nk,x0,2));
    otherwise
      error(['%s: %s, this type of driving function is not implemented ', ...
        'for a 2.5D plane wave.'],upper(mfilename),driving_functions);
  end
else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
