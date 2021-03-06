function status = test_secondary_source_positions(modus)
%TEST_SECONDARY_SOURCE_POSITIONS tests the correctness of the function
%secondary_source_positions()
%
%   Usage: status = test_secondary_source_positions(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%                 2: numerical verbose
%
%   Output parameters:
%       status - true or false
%
%   TEST_SECONDARY_SOURCE_POSITIONS(modus) checks if the function, that
%   calculates the secondary source positions and directions is working
%   correctly.

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


status = false;


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Main ============================================================
conf = SFS_config;
conf.secondary_sources.size = 4;
% reference values
x0_linear_ref = [
   -2.0000         0         0         0    -1.0000         0        0.2
   -1.8000         0         0         0    -1.0000         0        0.2
   -1.6000         0         0         0    -1.0000         0        0.2
   -1.4000         0         0         0    -1.0000         0        0.2
   -1.2000         0         0         0    -1.0000         0        0.2
   -1.0000         0         0         0    -1.0000         0        0.2
   -0.8000         0         0         0    -1.0000         0        0.2
   -0.6000         0         0         0    -1.0000         0        0.2
   -0.4000         0         0         0    -1.0000         0        0.2
   -0.2000         0         0         0    -1.0000         0        0.2
         0         0         0         0    -1.0000         0        0.2
    0.2000         0         0         0    -1.0000         0        0.2
    0.4000         0         0         0    -1.0000         0        0.2
    0.6000         0         0         0    -1.0000         0        0.2
    0.8000         0         0         0    -1.0000         0        0.2
    1.0000         0         0         0    -1.0000         0        0.2
    1.2000         0         0         0    -1.0000         0        0.2
    1.4000         0         0         0    -1.0000         0        0.2
    1.6000         0         0         0    -1.0000         0        0.2
    1.8000         0         0         0    -1.0000         0        0.2
    2.0000         0         0         0    -1.0000         0        0.2
    ];
x0_circle_ref = [
   2.00000   0.00000   0.00000  -1.00000   0.00000   0.00000   0.1995
   1.99006   0.19914   0.00000  -0.99503  -0.09957   0.00000   0.1995
   1.96034   0.39629   0.00000  -0.98017  -0.19815   0.00000   0.1995
   1.91115   0.58951   0.00000  -0.95557  -0.29476   0.00000   0.1995
   1.84295   0.77687   0.00000  -0.92148  -0.38843   0.00000   0.1995
   1.75644   0.95651   0.00000  -0.87822  -0.47825   0.00000   0.1995
   1.65248   1.12664   0.00000  -0.82624  -0.56332   0.00000   0.1995
   1.53209   1.28558   0.00000  -0.76604  -0.64279   0.00000   0.1995
   1.39647   1.43173   0.00000  -0.69824  -0.71587   0.00000   0.1995
   1.24698   1.56366   0.00000  -0.62349  -0.78183   0.00000   0.1995
   1.08509   1.68005   0.00000  -0.54255  -0.84003   0.00000   0.1995
   0.91242   1.77974   0.00000  -0.45621  -0.88987   0.00000   0.1995
   0.73068   1.86175   0.00000  -0.36534  -0.93087   0.00000   0.1995
   0.54168   1.92525   0.00000  -0.27084  -0.96262   0.00000   0.1995
   0.34730   1.96962   0.00000  -0.17365  -0.98481   0.00000   0.1995
   0.14946   1.99441   0.00000  -0.07473  -0.99720   0.00000   0.1995
  -0.04986   1.99938   0.00000   0.02493  -0.99969   0.00000   0.1995
  -0.24869   1.98448   0.00000   0.12434  -0.99224   0.00000   0.1995
  -0.44504   1.94986   0.00000   0.22252  -0.97493   0.00000   0.1995
  -0.63697   1.89585   0.00000   0.31849  -0.94793   0.00000   0.1995
  -0.82257   1.82301   0.00000   0.41129  -0.91151   0.00000   0.1995
  -1.00000   1.73205   0.00000   0.50000  -0.86603   0.00000   0.1995
  -1.16749   1.62388   0.00000   0.58374  -0.81194   0.00000   0.1995
  -1.32337   1.49956   0.00000   0.66169  -0.74978   0.00000   0.1995
  -1.46610   1.36035   0.00000   0.73305  -0.68017   0.00000   0.1995
  -1.59427   1.20761   0.00000   0.79713  -0.60380   0.00000   0.1995
  -1.70658   1.04287   0.00000   0.85329  -0.52144   0.00000   0.1995
  -1.80194   0.86777   0.00000   0.90097  -0.43388   0.00000   0.1995
  -1.87939   0.68404   0.00000   0.93969  -0.34202   0.00000   0.1995
  -1.93815   0.49351   0.00000   0.96908  -0.24676   0.00000   0.1995
  -1.97766   0.29808   0.00000   0.98883  -0.14904   0.00000   0.1995
  -1.99751   0.09969   0.00000   0.99876  -0.04985   0.00000   0.1995
  -1.99751  -0.09969   0.00000   0.99876   0.04985   0.00000   0.1995
  -1.97766  -0.29808   0.00000   0.98883   0.14904   0.00000   0.1995
  -1.93815  -0.49351   0.00000   0.96908   0.24676   0.00000   0.1995
  -1.87939  -0.68404   0.00000   0.93969   0.34202   0.00000   0.1995
  -1.80194  -0.86777   0.00000   0.90097   0.43388   0.00000   0.1995
  -1.70658  -1.04287   0.00000   0.85329   0.52144   0.00000   0.1995
  -1.59427  -1.20761   0.00000   0.79713   0.60380   0.00000   0.1995
  -1.46610  -1.36035   0.00000   0.73305   0.68017   0.00000   0.1995
  -1.32337  -1.49956   0.00000   0.66169   0.74978   0.00000   0.1995
  -1.16749  -1.62388   0.00000   0.58374   0.81194   0.00000   0.1995
  -1.00000  -1.73205   0.00000   0.50000   0.86603   0.00000   0.1995
  -0.82257  -1.82301   0.00000   0.41129   0.91151   0.00000   0.1995
  -0.63697  -1.89585   0.00000   0.31849   0.94793   0.00000   0.1995
  -0.44504  -1.94986   0.00000   0.22252   0.97493   0.00000   0.1995
  -0.24869  -1.98448   0.00000   0.12434   0.99224   0.00000   0.1995
  -0.04986  -1.99938   0.00000   0.02493   0.99969   0.00000   0.1995
   0.14946  -1.99441   0.00000  -0.07473   0.99720   0.00000   0.1995
   0.34730  -1.96962   0.00000  -0.17365   0.98481   0.00000   0.1995
   0.54168  -1.92525   0.00000  -0.27084   0.96262   0.00000   0.1995
   0.73068  -1.86175   0.00000  -0.36534   0.93087   0.00000   0.1995
   0.91242  -1.77974   0.00000  -0.45621   0.88987   0.00000   0.1995
   1.08509  -1.68005   0.00000  -0.54255   0.84003   0.00000   0.1995
   1.24698  -1.56366   0.00000  -0.62349   0.78183   0.00000   0.1995
   1.39647  -1.43173   0.00000  -0.69824   0.71587   0.00000   0.1995
   1.53209  -1.28558   0.00000  -0.76604   0.64279   0.00000   0.1995
   1.65248  -1.12664   0.00000  -0.82624   0.56332   0.00000   0.1995
   1.75644  -0.95651   0.00000  -0.87822   0.47825   0.00000   0.1995
   1.84295  -0.77687   0.00000  -0.92148   0.38843   0.00000   0.1995
   1.91115  -0.58951   0.00000  -0.95557   0.29476   0.00000   0.1995
   1.96034  -0.39629   0.00000  -0.98017   0.19815   0.00000   0.1995
   1.99006  -0.19914   0.00000  -0.99503   0.09957   0.00000   0.1995
   ];
x0_box_ref = [
    2.2000   -2.0000         0   -1.0000         0         0    0.2414
    2.2000   -1.8000         0   -1.0000         0         0    0.2000
    2.2000   -1.6000         0   -1.0000         0         0    0.2000
    2.2000   -1.4000         0   -1.0000         0         0    0.2000
    2.2000   -1.2000         0   -1.0000         0         0    0.2000
    2.2000   -1.0000         0   -1.0000         0         0    0.2000
    2.2000   -0.8000         0   -1.0000         0         0    0.2000
    2.2000   -0.6000         0   -1.0000         0         0    0.2000
    2.2000   -0.4000         0   -1.0000         0         0    0.2000
    2.2000   -0.2000         0   -1.0000         0         0    0.2000
    2.2000         0         0   -1.0000         0         0    0.2000
    2.2000    0.2000         0   -1.0000         0         0    0.2000
    2.2000    0.4000         0   -1.0000         0         0    0.2000
    2.2000    0.6000         0   -1.0000         0         0    0.2000
    2.2000    0.8000         0   -1.0000         0         0    0.2000
    2.2000    1.0000         0   -1.0000         0         0    0.2000
    2.2000    1.2000         0   -1.0000         0         0    0.2000
    2.2000    1.4000         0   -1.0000         0         0    0.2000
    2.2000    1.6000         0   -1.0000         0         0    0.2000
    2.2000    1.8000         0   -1.0000         0         0    0.2000
    2.2000    2.0000         0   -1.0000         0         0    0.2414
    2.0000    2.2000         0         0   -1.0000         0    0.2414
    1.8000    2.2000         0         0   -1.0000         0    0.2000
    1.6000    2.2000         0         0   -1.0000         0    0.2000
    1.4000    2.2000         0         0   -1.0000         0    0.2000
    1.2000    2.2000         0         0   -1.0000         0    0.2000
    1.0000    2.2000         0         0   -1.0000         0    0.2000
    0.8000    2.2000         0         0   -1.0000         0    0.2000
    0.6000    2.2000         0         0   -1.0000         0    0.2000
    0.4000    2.2000         0         0   -1.0000         0    0.2000
    0.2000    2.2000         0         0   -1.0000         0    0.2000
         0    2.2000         0         0   -1.0000         0    0.2000
   -0.2000    2.2000         0         0   -1.0000         0    0.2000
   -0.4000    2.2000         0         0   -1.0000         0    0.2000
   -0.6000    2.2000         0         0   -1.0000         0    0.2000
   -0.8000    2.2000         0         0   -1.0000         0    0.2000
   -1.0000    2.2000         0         0   -1.0000         0    0.2000
   -1.2000    2.2000         0         0   -1.0000         0    0.2000
   -1.4000    2.2000         0         0   -1.0000         0    0.2000
   -1.6000    2.2000         0         0   -1.0000         0    0.2000
   -1.8000    2.2000         0         0   -1.0000         0    0.2000
   -2.0000    2.2000         0         0   -1.0000         0    0.2414
   -2.2000    2.0000         0    1.0000         0         0    0.2414
   -2.2000    1.8000         0    1.0000         0         0    0.2000
   -2.2000    1.6000         0    1.0000         0         0    0.2000
   -2.2000    1.4000         0    1.0000         0         0    0.2000
   -2.2000    1.2000         0    1.0000         0         0    0.2000
   -2.2000    1.0000         0    1.0000         0         0    0.2000
   -2.2000    0.8000         0    1.0000         0         0    0.2000
   -2.2000    0.6000         0    1.0000         0         0    0.2000
   -2.2000    0.4000         0    1.0000         0         0    0.2000
   -2.2000    0.2000         0    1.0000         0         0    0.2000
   -2.2000         0         0    1.0000         0         0    0.2000
   -2.2000   -0.2000         0    1.0000         0         0    0.2000
   -2.2000   -0.4000         0    1.0000         0         0    0.2000
   -2.2000   -0.6000         0    1.0000         0         0    0.2000
   -2.2000   -0.8000         0    1.0000         0         0    0.2000
   -2.2000   -1.0000         0    1.0000         0         0    0.2000
   -2.2000   -1.2000         0    1.0000         0         0    0.2000
   -2.2000   -1.4000         0    1.0000         0         0    0.2000
   -2.2000   -1.6000         0    1.0000         0         0    0.2000
   -2.2000   -1.8000         0    1.0000         0         0    0.2000
   -2.2000   -2.0000         0    1.0000         0         0    0.2414
   -2.0000   -2.2000         0         0    1.0000         0    0.2414
   -1.8000   -2.2000         0         0    1.0000         0    0.2000
   -1.6000   -2.2000         0         0    1.0000         0    0.2000
   -1.4000   -2.2000         0         0    1.0000         0    0.2000
   -1.2000   -2.2000         0         0    1.0000         0    0.2000
   -1.0000   -2.2000         0         0    1.0000         0    0.2000
   -0.8000   -2.2000         0         0    1.0000         0    0.2000
   -0.6000   -2.2000         0         0    1.0000         0    0.2000
   -0.4000   -2.2000         0         0    1.0000         0    0.2000
   -0.2000   -2.2000         0         0    1.0000         0    0.2000
         0   -2.2000         0         0    1.0000         0    0.2000
    0.2000   -2.2000         0         0    1.0000         0    0.2000
    0.4000   -2.2000         0         0    1.0000         0    0.2000
    0.6000   -2.2000         0         0    1.0000         0    0.2000
    0.8000   -2.2000         0         0    1.0000         0    0.2000
    1.0000   -2.2000         0         0    1.0000         0    0.2000
    1.2000   -2.2000         0         0    1.0000         0    0.2000
    1.4000   -2.2000         0         0    1.0000         0    0.2000
    1.6000   -2.2000         0         0    1.0000         0    0.2000
    1.8000   -2.2000         0         0    1.0000         0    0.2000
    2.0000   -2.2000         0         0    1.0000         0    0.2414
    ];


%% ===== Test secondary source positions =================================
% Calculate current values
% linear array
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 21;
x0_linear = secondary_source_positions(conf);
% circular array
conf.secondary_sources.geometry = 'circle';
conf.secondary_sources.number = 63;
x0_circle = secondary_source_positions(conf);
% box form array
conf.secondary_sources.geometry = 'box';
conf.secondary_sources.number = 84;
x0_box = secondary_source_positions(conf);

if modus==0
    % Numerical mode (quiet)
    if ~all(eq(size(x0_linear),size(x0_linear_ref))) || ...
            ~all(eq(size(x0_circle),size(x0_circle_ref))) || ...
            ~all(eq(size(x0_box),size(x0_box_ref))) || ...
            ~all(abs(x0_linear(:)-x0_linear_ref(:))<1e-4) || ...
            ~all(abs(x0_circle(:)-x0_circle_ref(:))<1e-4) || ...
            ~all(abs(x0_box(:)-x0_box_ref(:))<1e-4)
        return;
    end
elseif modus==2
    if ~all(eq(size(x0_linear),size(x0_linear_ref)))
        error('%s: wrong size of linear array.',upper(mfilename));
    elseif ~all(abs(x0_linear(:)-x0_linear_ref(:))<1e-4)
        error('%s: wrong value at linear array.',upper(mfilename));
    end
    if ~all(eq(size(x0_circle),size(x0_circle_ref)))
        error('%s: wrong size of circular array.',upper(mfilename));
    elseif ~all(abs(x0_circle(:)-x0_circle_ref(:))<1e-4)
        error('%s: wrong value at circular array.',upper(mfilename));
    end
    if ~all(eq(size(x0_box),size(x0_box_ref))) 
        error('%s: wrong size of box shaped array.',upper(mfilename));
    elseif ~all(abs(x0_box(:)-x0_box_ref(:))<1e-4)
        error('%s: wrong value at box shaped array.',upper(mfilename));
    end
elseif modus==1
    % Graphical mode
    close all;
    % draw results
    figure
    title('Linear loudspeaker array');
    draw_loudspeakers(x0_linear,[1 1 0],conf);
    figure
    title('Circular loudspeaker array');
    draw_loudspeakers(x0_circle,[1 1 0],conf);
    figure
    title('Box shape loudspeaker array');
    draw_loudspeakers(x0_box,[1 1 0],conf);
else
    error(['%s: modus has to be 0 (numerical quiet), 1 (numerical), ', ...
            'or 2 (graphical).'],upper(mfilename));
end


status = true;
